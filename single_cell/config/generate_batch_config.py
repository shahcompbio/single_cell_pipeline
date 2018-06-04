
import os
import yaml
import json
import re
import single_cell
from single_cell.utils import helpers
import warnings


def get_version(reference):
    if reference.get("version", None):
        version = reference["version"].replace('.', '_')
    else:
        version = single_cell.__version__

        # strip setuptools metadata
        version = version.split("+")[0]

        # batch doesnt support . in pool names
        version = version.replace('.', '_')

    reference["version"] = version


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

    get_version(params)

    resolvedpaths = go_to_leafs(config, params)

    if output:
        helpers.makedirs(output, isfile=True)
        with open(output, 'w') as outfile:
            data = yaml.safe_dump(resolvedpaths, default_flow_style=False)
            outfile.write(data)
    else:
        return resolvedpaths


def generate_submit_config_in_temp(args):

    if args['which'] in ['clean_sentinels','generate_config']:
        return args

    if args.get("submit_config", None):
        return args

    batch_yaml = "batch.yaml"

    tmpdir = args.get("tmpdir", None)
    pipelinedir = args.get("pipelinedir", None)

    # use pypeliner tmpdir to store yaml
    if pipelinedir:
        batch_yaml = os.path.join(pipelinedir, batch_yaml)
    elif tmpdir:
        batch_yaml = os.path.join(tmpdir, batch_yaml)
    else:
        warnings.warn("no tmpdir specified, generating configs in working dir")
        batch_yaml = os.path.join(os.getcwd(), batch_yaml)

    batch_yaml = helpers.get_incrementing_filename(batch_yaml)

    params_override = args["config_override"]

    main(output=batch_yaml, input_params=params_override)

    args["submit_config"] = batch_yaml

    return args

if __name__ == "__main__":
    input_params = '{"cluster": "azure", "aligner": "bwa-aln", "reference": "grch37", "smoothing_function": "modal"}'
    input_params = json.loads(input_params)

    main(input_params=input_params, output="out.yaml")
