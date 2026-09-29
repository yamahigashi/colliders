"""Capture the vertical first-pass surface using the HEAD plugin."""

import json

from helpers import maya_session, parser


def main():
    arguments = parser(__doc__)
    arguments.add_argument(
        "--build-configuration",
        required=True,
        help="Compiler, target and configuration shared by the reference and candidate builds",
    )
    args = arguments.parse_args()
    with maya_session(args.plugin) as metadata:
        import maya.cmds as cmds
        from test_continuous_lift_activation import LiftFixture, REST_REFERENCE, packed_points

        fixture = LiftFixture()
        fixture.pose(0.0)
        fixture.select_rows(fixture.base(), (10.0,))
        dimensions, points = fixture.read()
        providers = [
            plugin
            for plugin in cmds.pluginInfo(query=True, listPlugins=True)
            if "yddSkirtBellCollider" in (cmds.pluginInfo(plugin, query=True, dependNode=True) or [])
        ]
        if len(providers) != 1:
            raise RuntimeError("Expected one plugin providing yddSkirtBellCollider")
        reference = dict(
            dimensions=list(dimensions),
            cv_xyzw_le_hex=packed_points(points).hex(),
            maya_version=metadata["maya"],
            maya_api_version=cmds.about(apiVersion=True),
            plugin_version=cmds.pluginInfo(providers[0], query=True, version=True),
            plugin_sha256=metadata["sha256"],
            python_version=metadata["python"],
            build_configuration=args.build_configuration,
            capture="Normal outputSurface at theta zero in a fresh synthetic fixture",
        )
        REST_REFERENCE.parent.mkdir(parents=True, exist_ok=True)
        REST_REFERENCE.write_text(json.dumps(reference, indent=2) + "\n", encoding="utf-8")
        print("Wrote", REST_REFERENCE.name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
