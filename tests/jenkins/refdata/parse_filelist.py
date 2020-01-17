import json

with open('files.json', "r") as input, open("files.txt", "w") as output:
    data = json.load(input)
    for f in data:
        output.write("%s\n" % f["name"])