import C17_geometry_arcs as c17

# run file
jsonfile = 'input_c17_arcs.json'
info = c17.read_json(jsonfile)
for i in range(len(info["run"])):
    c17.main(info,info["run"][i])