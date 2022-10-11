import machupX as mx

wing = mx.Scene("wingscene.json")

# wing.display_wireframe(show_vortices=False)

wing.export_dxf(aircraft="wing",number_guide_curves=6,export_as_prismoid=True)

