import C17_geometry as c17

# run file
jsonfile = 'input_c17.json'
info = c17.read_json(jsonfile)
for i in range(len(info["run"])):
    c17.main(info,info["run"][i])