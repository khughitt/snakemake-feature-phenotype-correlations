"""
Feature-phenotype correlation pipeline
KH Feb 2020
"""
import os

out_dir = os.path.join(config["output_dir"], config["name"], config["version"])

wildcard_constraints:
    feature="({})".format(")|(".join(config["features"].keys())),
    phenotype="({})".format(")|(".join(config["phenotypes"].keys())),
    method="({})".format(")|(".join(config["methods"]))

rule all:
    input:
        expand(os.path.join(out_dir, "correlations/features/{feature}_{method}.feather"),
               feature=config["features"].keys(), method=config["methods"]),
        expand(os.path.join(out_dir, "correlations/phenotypes/{phenotype}_{method}.feather"),
               phenotype=config["phenotypes"].keys(), method=config["methods"]),
        expand(os.path.join(out_dir, "correlations/feature-phenotype/{feature}_{phenotype}_{method}.feather"),
               feature=config["features"].keys(), phenotype=config["phenotypes"].keys(),
               method=config["methods"])

rule compute_cross_correlations:
    input:
        features=os.path.join(out_dir, "input/features/{feature}.feather"),
        phenotypes=os.path.join(out_dir, "input/phenotypes/{phenotype}.feather")
    output:
        os.path.join(out_dir, "correlations/feature-phenotype/{feature}_{phenotype}_{method}.feather")
    script:
        "src/compute_cross_correlations.py"

rule compute_phenotype_correlations:
    input:
        os.path.join(out_dir, "input/phenotypes/{phenotype}.feather")
    output:
        os.path.join(out_dir, "correlations/phenotypes/{phenotype}_{method}.feather")
    script:
        "src/compute_correlations.py"

rule compute_feature_correlations:
    input:
        os.path.join(out_dir, "input/features/{feature}.feather")
    output:
        os.path.join(out_dir, "correlations/features/{feature}_{method}.feather")
    script:
        "src/compute_correlations.py"

rule create_phenotype_symlinks:
    output:
        os.path.join(out_dir, "input/phenotypes/{phenotype}.feather")
    run:
        os.symlink(config["phenotypes"][wildcards["phenotype"]], output[0])

rule create_feature_symlinks:
    output:
        os.path.join(out_dir, "input/features/{feature}.feather")
    run:
        os.symlink(config["features"][wildcards["feature"]], output[0])
