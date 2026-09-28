#pragma once

#include <array>
#include <cstddef>
#include <list>
#include <memory>
#include <mutex>

// Payload counts the objects and their vector capacities. Borrowed pointers can
// outlive eviction; those objects are no longer part of the node's retained cache.
template<class Topology, class Preparation>
class SkirtSurfaceCache
{
public:
    using TopologyPtr = std::shared_ptr<const Topology>;
    using PreparationPtr = std::shared_ptr<const Preparation>;
    static const std::size_t maxTopologies = 4;
    static const std::size_t maxPreparations = 2;
    static const std::size_t maxPayloadBytes = 4 * 1024 * 1024;

    template<class Equal>
    TopologyPtr findTopology(const Topology& key, Equal equal)
    {
        for (;;)
        {
            const auto snapshot = topologySnapshot();
            TopologyPtr found;
            for (std::size_t i = 0; i < snapshot.count; ++i)
                if (equal(*snapshot.values[i], key))
                {
                    found = snapshot.values[i];
                    break;
                }
            std::lock_guard<std::mutex> lock(mutex_);
            if (revision_ != snapshot.revision)
                continue;
            if (found)
            {
                auto it = findEntry(found);
                entries_.splice(entries_.begin(), entries_, it);
                ++revision_;
            }
            return found;
        }
    }

    template<class Equal>
    TopologyPtr publishTopology(const TopologyPtr& built, std::size_t payloadBytes, Equal equal)
    {
        // Allocate the list node before taking the node cache lock.
        std::list<Entry> staged;
        if (payloadBytes <= maxPayloadBytes)
            staged.emplace_back(built, payloadBytes);
        std::list<Entry> retired;
        for (;;)
        {
            const auto snapshot = topologySnapshot();
            TopologyPtr found;
            for (std::size_t i = 0; i < snapshot.count; ++i)
                if (equal(*snapshot.values[i], *built))
                {
                    found = snapshot.values[i];
                    break;
                }
            {
                std::lock_guard<std::mutex> lock(mutex_);
                if (revision_ != snapshot.revision)
                    continue;
                if (found)
                {
                    auto it = findEntry(found);
                    entries_.splice(entries_.begin(), entries_, it);
                    ++revision_;
                    return found;
                }
                if (!staged.empty())
                {
                    entries_.splice(entries_.begin(), staged, staged.begin());
                    retainedBytes_ += payloadBytes;
                    while (entries_.size() > maxTopologies || retainedBytes_ > maxPayloadBytes)
                    {
                        auto victim = --entries_.end();
                        retainedBytes_ -= victim->payloadBytes;
                        retired.splice(retired.end(), entries_, victim);
                    }
                    ++revision_;
                }
            }
            return built;
        }
    }

    template<class Matches>
    PreparationPtr findPreparation(const TopologyPtr& topology, Matches matches)
    {
        for (;;)
        {
            const auto snapshot = preparationSnapshot(topology);
            if (!snapshot.entryPresent)
                return {};
            PreparationPtr found;
            for (std::size_t i = 0; i < snapshot.count; ++i)
                if (matches(*snapshot.values[i]))
                {
                    found = snapshot.values[i];
                    break;
                }
            std::lock_guard<std::mutex> lock(mutex_);
            if (revision_ != snapshot.revision)
                continue;
            if (found)
            {
                auto entry = findEntry(topology);
                auto prep = findPreparationEntry(entry, found);
                entry->preparations.splice(entry->preparations.begin(), entry->preparations, prep);
                entries_.splice(entries_.begin(), entries_, entry);
                ++revision_;
            }
            return found;
        }
    }

    template<class Matches>
    PreparationPtr publishPreparation(const TopologyPtr& topology, const PreparationPtr& built,
                                      std::size_t payloadBytes, Matches matches)
    {
        std::list<PreparationEntry> staged;
        if (payloadBytes <= maxPayloadBytes)
            staged.emplace_back(built, payloadBytes);
        std::list<PreparationEntry> retiredPreparations;
        std::list<Entry> retiredTopologies;
        for (;;)
        {
            const auto snapshot = preparationSnapshot(topology);
            if (!snapshot.entryPresent)
                return built;
            PreparationPtr found;
            for (std::size_t i = 0; i < snapshot.count; ++i)
                if (matches(*snapshot.values[i]))
                {
                    found = snapshot.values[i];
                    break;
                }
            {
                std::lock_guard<std::mutex> lock(mutex_);
                if (revision_ != snapshot.revision)
                    continue;
                auto entry = findEntry(topology);
                if (found)
                {
                    auto prep = findPreparationEntry(entry, found);
                    entry->preparations.splice(entry->preparations.begin(), entry->preparations, prep);
                    entries_.splice(entries_.begin(), entries_, entry);
                    ++revision_;
                    return found;
                }
                if (staged.empty() || entry->payloadBytes + payloadBytes > maxPayloadBytes)
                    return built;
                entry->preparations.splice(entry->preparations.begin(), staged, staged.begin());
                retainedBytes_ += payloadBytes;
                entry->payloadBytes += payloadBytes;
                entries_.splice(entries_.begin(), entries_, entry);
                while (entries_.front().preparations.size() > maxPreparations)
                    retirePreparation(entries_.begin(), retiredPreparations);
                while (retainedBytes_ > maxPayloadBytes && entries_.size() > 1)
                {
                    auto victim = --entries_.end();
                    retainedBytes_ -= victim->payloadBytes;
                    retiredTopologies.splice(retiredTopologies.end(), entries_, victim);
                }
                while (retainedBytes_ > maxPayloadBytes && entries_.front().preparations.size() > 1)
                    retirePreparation(entries_.begin(), retiredPreparations);
                ++revision_;
            }
            return built;
        }
    }

    std::size_t retainedPayloadBytes() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return retainedBytes_;
    }

    std::size_t topologyCount() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        return entries_.size();
    }

    std::size_t preparationCount() const
    {
        std::lock_guard<std::mutex> lock(mutex_);
        std::size_t count = 0;
        for (const auto& entry : entries_)
            count += entry.preparations.size();
        return count;
    }

private:
    struct PreparationEntry
    {
        PreparationEntry(const PreparationPtr& value, std::size_t bytes) : value(value), bytes(bytes) {}
        PreparationPtr value;
        std::size_t bytes;
    };
    struct Entry
    {
        Entry(const TopologyPtr& value, std::size_t bytes) : value(value), payloadBytes(bytes) {}
        TopologyPtr value;
        std::size_t payloadBytes;
        std::list<PreparationEntry> preparations;
    };
    struct TopologySnapshot
    {
        std::array<TopologyPtr, maxTopologies> values;
        std::size_t count = 0;
        std::size_t revision = 0;
    };
    struct PreparationSnapshot
    {
        std::array<PreparationPtr, maxPreparations> values;
        std::size_t count = 0;
        std::size_t revision = 0;
        bool entryPresent = false;
    };

    TopologySnapshot topologySnapshot() const
    {
        TopologySnapshot snapshot;
        std::lock_guard<std::mutex> lock(mutex_);
        snapshot.revision = revision_;
        for (const auto& entry : entries_)
            snapshot.values[snapshot.count++] = entry.value;
        return snapshot;
    }

    PreparationSnapshot preparationSnapshot(const TopologyPtr& topology) const
    {
        PreparationSnapshot snapshot;
        std::lock_guard<std::mutex> lock(mutex_);
        snapshot.revision = revision_;
        for (const auto& entry : entries_)
            if (entry.value == topology)
            {
                snapshot.entryPresent = true;
                for (const auto& preparation : entry.preparations)
                    snapshot.values[snapshot.count++] = preparation.value;
                break;
            }
        return snapshot;
    }

    typename std::list<Entry>::iterator findEntry(const TopologyPtr& topology)
    {
        for (auto it = entries_.begin(); it != entries_.end(); ++it)
            if (it->value == topology)
                return it;
        return entries_.end();
    }

    typename std::list<PreparationEntry>::iterator findPreparationEntry(
        typename std::list<Entry>::iterator entry, const PreparationPtr& preparation)
    {
        for (auto it = entry->preparations.begin(); it != entry->preparations.end(); ++it)
            if (it->value == preparation)
                return it;
        return entry->preparations.end();
    }

    void retirePreparation(typename std::list<Entry>::iterator entry, std::list<PreparationEntry>& retired)
    {
        auto victim = --entry->preparations.end();
        retainedBytes_ -= victim->bytes;
        entry->payloadBytes -= victim->bytes;
        retired.splice(retired.end(), entry->preparations, victim);
    }

    mutable std::mutex mutex_;
    std::list<Entry> entries_;
    std::size_t retainedBytes_ = 0;
    std::size_t revision_ = 0;
};
