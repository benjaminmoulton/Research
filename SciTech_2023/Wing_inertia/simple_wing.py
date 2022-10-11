import machupX as mx

wings = mx.Scene("wingscene.json")

# wings.display_wireframe(show_vortices=False)

wings.export_dxf(aircraft="wing",number_guide_curves=6,export_as_prismoid=False)

