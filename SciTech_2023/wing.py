import machupX as mx

twist = mx.Scene("twistscene.json")

# twist.display_wireframe(show_vortices=False)

twist.export_dxf(aircraft="twist",number_guide_curves=6,export_as_prismoid=True)

