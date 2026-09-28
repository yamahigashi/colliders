# ADR-0003: ydd プラグインの識別子

Status: Accepted（2026-09-10 実装・配置済み）
Date: 2026-09-10
Owner: yamahigashi colliders maintainers
Supersedes: 既存 ADR のプラグイン名と Maya ノード型名、ノード ID の表記
Superseded by: none

## 背景

fork 元と本 fork は `colliders.mll`、3種類のノード名と ID、描画登録、Python モジュール、AE テンプレートを共有している。
利用者は重複の解消と `ydd` 接頭辞を指定し、割り当て済みの `MTypeId` 範囲は持っていないと回答した。
ノード名だけを変更しても ID とプラグイン名が重複するため、Maya が認識する識別子をまとめて分離する。

## ノード名と ID

本 fork の5ノードを以下に固定する。
番号は内部利用向け範囲 `0x00000000..0x0007ffff` 内の、このリポジトリで管理するローカル割当である。
Autodesk が割り当てたグローバル ID ではない。

| 旧ノード型 | 新ノード型 | 旧 ID | 新 ID |
| --- | --- | ---: | --- |
| `bellCollider` | `yddBellCollider` | 1274434 | `0x0007DD00` |
| `planeCollider` | `yddPlaneCollider` | 1274435 | `0x0007DD01` |
| `skirtBellCollider` | `yddSkirtBellCollider` | 1274436 | `0x0007DD02` |
| `skirtCollideDeformer` | `yddSkirtCollideDeformer` | 1274437 | `0x0007DD03` |
| `skirtWaveDeformer` | `yddSkirtWaveDeformer` | 1274438 | `0x0007DD04` |

2026-09-15 追記：面フィット用の依存ノードを 1 個追加し、同じローカル範囲から採番する。

| ノード型 | 種別 | ID |
| --- | --- | --- |
| `yddSkirtSurfaceFit` | `MPxNode`（`kDependNode`、描画分類なし） | `0x0007DD05` |

名前と ID の C++ 定義を `sources/pluginIdentity.h` に集約し、登録処理と各クラスで参照する。
実行時の採番や環境変数による ID の切り替えは行わない。
手元のソースを検索し、Maya 内でも fork 元との同時登録を検証する。
この検証は別環境のプラグインとの非衝突を保証しない。

Autodesk は内部利用向け範囲と配布向けのグローバル割当を区別している。
組織をまたぐ配布へ移行する場合は、登録済み範囲の取得と保存済みシーンの移行を一緒に設計する。
[MTypeId の公式説明](https://help.autodesk.com/cloudhelp/2026/ENU/MAYA-API-REF/cpp_ref/class_m_type_id.html)

## プラグインとスクリプトの識別子

ビルドターゲットと Maya のプラグイン名は `yddColliders`、Windows バイナリは `yddColliders.mll` とする。
プラグインのバージョンは `3.0.0` とする。
リポジトリ名と既存のパッケージディレクトリは変更しない。

Python の公開モジュールは `scripts/yddColliders.py` とし、利用側を `import yddColliders` に変更する。
公開関数 `createBellCollider`、`createSkirtBellCollider`、`show` の役割と引数は維持する。
UI のウィンドウ ID と、これらの関数が作るカスタムノードの既定名にも `ydd` を付ける。

ロケーターの描画分類はそれぞれ `drawdb/geometry/yddBellCollider`、`drawdb/geometry/yddPlaneCollider`、`drawdb/geometry/yddSkirtBellCollider` とする。
描画登録 ID は `yddCollidersPlugin` とし、登録と解除で同じ値を使う。
AE テンプレートはファイル名と global proc 名を `AEyddBellColliderTemplate`、`AEyddSkirtBellColliderTemplate` に揃える。

旧名のノード登録、旧 ID の登録、旧 Python モジュールの別名は残さない。
fork 元は旧名を継続して使用できる。
C++ の実装クラス名と既存のソースファイル名は内部識別子なので維持する。

## 変形と利用側の契約

ノードの属性名、型、既定値、変形計算、スケジューリングは変更しない。
ノード名と ID の分離前後で、同じ入力の形状出力を完全一致させる。
既存の衝突、リング、Wave の制約も継続する。

利用側のスカートコンポーネントは `yddColliders` をロードし、`ydd` を付けたノードだけを作成する。
fork 元が先にロードされていても利用するプラグインを変えない。
Component と Guide のバージョンはともに `[3, 0, 0]` とする。
ガイドのパラメータとホスト属性は変更しない。

## 保存済みシーンと配置

型名と ID の変更は、既存シーンの互換性を切る変更として扱う。
ID は `.mb` に保存されるため、旧シーンを新プラグインだけでそのまま読むことは保証しない。
旧シーンは旧バイナリで開き、ガイドから新しいリグを構築するか、必要なデータを別途移行する。
この変更で利用者のシーンを自動変換しない。
[MTypeId とバイナリ保存の契約](https://help.autodesk.com/cloudhelp/2026/ENU/MAYA-API-REF/cpp_ref/class_m_type_id.html)

変更前のバイナリと検証コードを `build/identity-validation/baseline/` に保存する。
新バイナリの配置先は `plugins/<Maya>/yddColliders.mll` と、コンポーネント側の `platforms/win64_<Maya>/plug-ins/yddColliders.mll` とする。
検証後に、両配置先にある旧 fork の `colliders.mll` を検索パスから除く。
fork 元のリポジトリと、コンポーネント側の `platforms_new`、`platforms_old` など利用者の別保管場所は変更しない。

過去の ADR、検証結果、画像採取メタデータは当時の識別子を残す。
現在の README、実行コード、テスト、ビルド手順は新しい識別子に更新する。

## 検討した選択肢

| 選択肢 | 判断 |
| --- | --- |
| ノード名だけを変更する | ID とプラグイン名の重複が残るため採用しない |
| 旧 ID の近くの番号へ変更する | その範囲の所有と空きを確認できないため採用しない |
| グローバル ID の取得まで実装を止める | 今回は内部利用向けの固定 ID で進め、配布時の条件を明示する |
| 旧名も同時登録する | fork 元との衝突を再導入するため採用しない |
| 内部利用 ID と ydd 識別子に統一する | 現在の重複を解消できるため採用する |

## 検収条件

- Maya 2022 から2027の各 SDK でビルドし、各版で既存テストと名前、ID、描画分類の検証を実行する。
- 変更前の fork と形状出力を比較する。
- fork 元を先にロードする順序と、本 fork を先にロードする順序で、両方のノードを作成できることを確認する。
- `.ma` と `.mb` を保存して開き直し、新しいノード型と ID、属性、接続が維持されることを確認する。
- 本 fork の解除が fork 元のノード登録に影響しないことを確認する。
- 実際の mGear リグを構築し、新しいプラグインとノードの利用、Wave の既定値と出力を確認する。
- 配置した各バイナリのハッシュを記録し、通常の配置先から旧 fork のバイナリを除いたことを確認する。

## 再検討の条件

外部配布、別の内部プラグインとの ID 衝突、既存シーンの一括移行が必要になった場合に、この識別子と移行の契約を見直す。

## 実装と検収

2026-09-10に実装と配置を完了した。
Maya 2022 から2027での共存、保存後の再読込、変形の一致、実リグの検証結果は [検収記録](../validation/ydd-identity.md) に記載した。
