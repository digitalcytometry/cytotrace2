import argparse


def argument_parser():
    parser = argparse.ArgumentParser(description="CytoTRACE 2 is a computational method for "
                                                 "predicting cellular potency categories and "
                                                 "absolute developmental potential from "
                                                 "single-cell RNA-sequencing data.")

    # Required arguments
    required = parser.add_argument_group('Required arguments')
    required.add_argument("-f", "--input-path", help="Path to input scRNA-seq data", type=str,
                          default=None, required=True)


    # Options

    parser.add_argument("-a", "--annotation-path", type=str, default="",
                        help="Path to input annotations")

    parser.add_argument("-sp", "--species", default="mouse",
                        help="Species name from which the gene names for the input data come, default 'mouse'",
                        choices=["mouse", "human"])

    parser.add_argument("-fm", "--full-model", help="Use full ensemble of 17 models rather than default reduced ensemble of 5 most predictive models", action="store_true")

    parser.add_argument("-bs", "--batch-size", help="Integer or NULL, indicating the number of cells to subsample for the pipeline steps. No subsampling if NULL. (default is 10000, recommended value for input data size > 10K cells)", type=int, default=10000)    

    parser.add_argument("-sbs", "--smooth-batch-size", help="Integer or NULL, indicating the number of cells to subsample further for the smoothing step of the pipeline. No subsampling if NULL. (default is 1000, recommended value for input data size > 1K cells)", type=int, default=1000)    

    parser.add_argument("-dp", "--disable-parallelization", help="Disable parallelization for CytoTRACE 2", action="store_true")

    parser.add_argument("-mc", "--max-cores", help="Integer, indicating the maximum number of CPU cores to use for parallelization", type=int, default=None)  

    parser.add_argument("-mpc", "--max-pcs", help="Integer, indicating the maximum number of principal components to use in the smoothing by kNN step (default is 200)", type=int, default=200)  

    parser.add_argument("-r", "--seed", help="Integer, specifying the seed for reproducibility in random processes (default is 14).", type=int, default=14)  

    parser.add_argument("-o", "--output-dir", type=str, default="cytotrace2_results",
                        help="Directory path for output predictions and intermediate files (default is 'cytotrace2_results')")




    arguments = parser.parse_args()

    return arguments.__dict__
