
import os
import yaml
import json
import re
import copy


def get_version(reference):
    if reference.get("version", None):
        return reference["version"]

    version = single_cell.__version__

    # strip setuptools metadata
    version = version.split("+")[0]

    # batch doesnt support . in pool names
    version = version.replace('.', '_')

    return version


def eval_expression(fields, reference):

    if not isinstance(fields, str):
        return fields

    for match in re.finditer("\$\{(.*?)\}", fields):

        match = match.group()

        assert match.startswith("${") and match.endswith("}")

        reference_key = match[2:-1]

#         print reference_key, match, reference[reference_key]

        assert reference_key in reference

        fields = fields.replace(match, reference[reference_key])

    return fields


def go_to_leafs(datadict, reference):
    def recurse_dict(data, reference):
        '''
        recurse into dict and yield ylims.
        '''
        for key, val in data.iteritems():
            del data[key]
            newkey = eval_expression(key, reference)
            data[newkey] = val

            if isinstance(val, dict):
                recurse_dict(val, reference)
            else:
                if isinstance(val, list):
                    data[key] = [
                        eval_expression(
                            element,
                            reference) for element in val]
                elif isinstance(val, basestring):
                    data[key] = eval_expression(val, reference)
                elif type(val) in [int, float, bool] or not val:
                    continue
                else:
                    print val
                    print type(val)
                    raise Exception("Unknown {}".format(val))

    recurse_dict(datadict, reference)

    return datadict


def main(output=None, input_params=None):

    cfgdir = os.path.realpath(os.path.dirname(__file__))
    config = os.path.join(cfgdir, "batch.yaml")
    params = os.path.join(cfgdir, "batch_params.yaml")

    config = yaml.load(open(config))

    params = yaml.load(open(params))

    if input_params:
        params.update(input_params)

    resolvedpaths = go_to_leafs(config, params)

    if output:
        with open(output, 'w') as outfile:
            data = yaml.safe_dump(resolvedpaths, default_flow_style=False)
            outfile.write(data)
    else:
        return resolvedpaths

if __name__ == "__main__":
    input_params = '{"cluster": "azure", "aligner": "bwa-aln", "reference": "grch37", "smoothing_function": "modal"}'
    input_params = json.loads(input_params)

    main(input_params=input_params, output="out.yaml")
