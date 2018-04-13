import yaml
import json
import copy
import os
import single_cell


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

    if not fields.startswith('{') or not fields.endswith('}'):
        return fields

    # remove braces
    fields = fields[1:-1]
    fields = fields.split(":")

    orig_ref = copy.deepcopy(reference)

    for field in fields:
        field = orig_ref[field[1:]] if field.startswith('$') else field

        if isinstance(reference, dict):
            reference = reference[field]
        else:
            raise Exception("Could not resolve: {}".format(':'.join(fields)))

    return reference


def go_to_leafs(datadict, reference):
    def recurse_dict(data, reference):
        '''
        recurse into dict and yield ylims.
        '''
        for key, val in data.iteritems():
            if isinstance(val, dict):
                recurse_dict(val, reference)
            else:
                if isinstance(val, list):
                    data[key] = [
                        eval_expression(
                            element,
                            reference) for element in val]
                elif isinstance(val, str):
                    data[key] = eval_expression(val, reference)
                elif type(val) in [int, float, bool] or not val:
                    continue
                else:
                    raise Exception("Unknown {}".format(val))

    recurse_dict(datadict, reference)

    return datadict


def add_version_to_pools(cfgdict, reference):

    version = get_version(reference)

    if not cfgdict.get("pools", None):
        return cfgdict

    for poolid, poolname in cfgdict["pools"].iteritems():
        cfgdict["pools"][poolid] = '_'.join((poolname, version))

    return cfgdict


def main(output=None, input_params=None):

    cfgdir = os.path.realpath(os.path.dirname(__file__))
    config = os.path.join(cfgdir, "config.yaml")
    params = os.path.join(cfgdir, "config_params.yaml")

    config = yaml.load(open(config))

    params = yaml.load(open(params))

    if input_params:
        params.update(input_params)

    resolvedpaths = go_to_leafs(config, params)

    resolvedpaths = add_version_to_pools(resolvedpaths, params)

    if output:
        with open(output, 'w') as outfile:
            data = yaml.dump(resolvedpaths, default_flow_style=False)
            outfile.write(data)
    else:
        return resolvedpaths

if __name__ == "__main__":

    # DEBUG ONLY
    input_params = '{"cluster": "azure", "aligner": "bwa-aln", "reference": "grch37", "smoothing_function": "modal"}'
    input_params = json.loads(input_params)

    config = "config.yaml"
    params = "params.yaml"

    print main(input_params=input_params)
