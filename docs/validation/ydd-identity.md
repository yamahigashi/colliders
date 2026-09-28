# ydd 識別子の検収記録

検証日: 2026-09-10

[ADR-0003](../adr/0003-ydd-plugin-identity.md) に従い、プラグイン、5ノードの型名と ID、描画登録、Python モジュール、AE テンプレートを分離した。
`yddColliders.mll` はバージョン `3.0.0` で、利用側のスカートコンポーネントも同じ版で検証した。
名前と ID は [pluginIdentity.h](../../sources/pluginIdentity.h) に集約した。

## Maya 各版の検証

Windows x64 の Maya 2022 から2027について、それぞれの SDK で Release ビルドし、対応する `mayapy` で既存テストと登録検証を実行した。
結果、使用バイナリ、検証コードの SHA-256 はリポジトリ外の検証データに保存した。

| Maya | 実行数 | 成功 | Skip | 共存のロード順 | 保存後の再読込 |
| --- | ---: | ---: | ---: | --- | --- |
| 2022 | 28 | 27 | 1 | 両順序で成功 | `.ma`、`.mb` とも成功 |
| 2023 | 28 | 27 | 1 | 両順序で成功 | `.ma`、`.mb` とも成功 |
| 2024 | 28 | 27 | 1 | 両順序で成功 | `.ma`、`.mb` とも成功 |
| 2025 | 28 | 27 | 1 | 両順序で成功 | `.ma`、`.mb` とも成功 |
| 2026 | 28 | 28 | 0 | 両順序で成功 | `.ma`、`.mb` とも成功 |
| 2027 | 28 | 27 | 1 | 両順序で成功 | `.ma`、`.mb` とも成功 |

Skip は Maya 2026 専用の既存の座標ハッシュ比較である。
有限値、点数、周期面のチェックは全版で実行した。
追加した3テストは、5ノードの型名、ID、API 種別、プラグインの vendor と version、描画分類、Python モジュールと生成ノード名を検証する。

共存検証では、fork 元を先にロードする順序と本 fork を先にロードする順序を、別プロセスで実行した。
各回で計8ノードを作成し、8接続、属性、ID、出力を `.ma` と `.mb` の保存前後で比較した。
12回の共存検証と24回の再読込が成功し、出力差の最大値はすべて0だった。
一方のプラグインを強制指定なしで解除した後も他方を利用でき、Python モジュールの同時 import と AE テンプレートの MEL 手続き解決も全版で成功した。

fork 元の既存バイナリは、Maya 2022 と2023では2ノードのみを登録し、2024版は存在しなかった。
この3版では、現在の fork 元のソースを変更せずに別のビルド出力先でコンパイルし、3ノードを持つバイナリで検証した。
2025から2027では fork 元の既存バイナリを使った。

## 変形と実リグの一致

Maya 2026で、変更直前の fork のコードとバイナリを使った出力採取と、新しい識別子での出力採取を比較した。
275条件すべてで double 座標が完全一致し、最大絶対差は0だった。
比較対象には、衝突、Wave の信号、ゼロ寄与、DG、Serial、Parallel、時間の往復が含まれる。

配置済みバイナリでは、mGear 5.3.2を使ってガイドから3構成の実リグを構築した。
構成は Wave 無効で post-collision 有効、Wave 有効で post-collision 無効、両方有効である。
各リグの FK は40個、リングは2個で、新しいノード型、7ホスト属性の契約と直接接続、変形順、Wave の初期ゼロ変形、操作後の変形を確認した。
識別子変更前に採取した実リグの8出力ハッシュとも一致した。

利用側の自動ロードも、fork 元だけがロードされた状態から確認した。
実際の Component のプラグイン確認処理が `yddColliders` を別にロードし、両者の登録とロードパスを確認できた。
3構成とも両プラグインをロードしたまま構築した。

## 配置と静的検査

各版のバイナリを `plugins/<Maya>/yddColliders.mll` とコンポーネントの `platforms/win64_<Maya>/plug-ins/yddColliders.mll` に配置し、計12ファイルが検証済みビルドと同じ SHA-256 であることを確認した。
これらの配置先にあった旧 fork の `colliders.mll` は、ハッシュを確認して保存した後に除去した。
変更前のコードとバイナリは `build/identity-validation/baseline/` に保存し、55ファイルのハッシュを再確認した。

Ruff の lint は実行用 Python スクリプト、テスト、Component、Guide で成功し、テスト用15ファイルの format check と Component、Guide の ty 検査も成功した。
名前を変更した公開スクリプトと Component の既存文字列2か所には変更前から整形差分があるため、一括整形は行っていない。
両リポジトリの `git diff --check` は、Windows バッチの CRLF を許容する `cr-at-eol` 設定で成功した。

## 再現手順と範囲

ビルドは `buildAll.bat <Maya>`、テストは対応する `mayapy` で実行する。
共存検証には同じ Maya 版に対応する fork 元のバイナリとスクリプトを指定する。

```powershell
$mayapy = 'C:\Program Files\Autodesk\Maya2026\bin\mayapy.exe'
& $mayapy -B tests/run_maya_tests.py --plugin C:\fork\plugins\2026\yddColliders.mll
& $mayapy -B tests/identity/coexistence.py --plugin C:\fork\plugins\2026\yddColliders.mll --upstream-plugin C:\upstream\plugins\2026\colliders.mll --upstream-scripts C:\upstream\scripts --load-order upstream-first --output coexist-upstream-first.json
& $mayapy -B tests/identity/coexistence.py --plugin C:\fork\plugins\2026\yddColliders.mll --upstream-plugin C:\upstream\plugins\2026\colliders.mll --upstream-scripts C:\upstream\scripts --load-order fork-first --output coexist-fork-first.json
```

変更前との出力比較と実リグ検証の引数は [テスト手順](../../tests/README.md) に記載した。
生ログと座標データはリポジトリ外の `build/identity-validation/` にある。
実リグと変更前後の座標比較は Maya 2026で実行した。
今回の識別子変更では GUI の画素比較と Cached Playback の生成、復元は再検証していない。

型名と ID の変更により、旧シーンは旧バイナリで開き、ガイドからの再構築または別途データ移行が必要になる。
新しい ID `0x0007DD00..0x0007DD04` は内部利用向けのローカル割当で、Autodesk のグローバル割当ではない。
今回検証した fork 元以外のプラグインとの非衝突は保証しない。
ID の割当範囲とバイナリ保存の仕様は [Autodesk の MTypeId 説明](https://help.autodesk.com/cloudhelp/2026/ENU/MAYA-API-REF/cpp_ref/class_m_type_id.html) に従う。
