library(optparse)

parseArgs <- function(){
  # Parses in command line arguments
  #
  # Args:
  #   NA
  #
  # Returns:
  #   Dataframe with the columns ntrees filled with
  #   value supplied by the command line
  #Load in the arguments from the command line
  option_list = list(
    make_option(c("-g", "--growthrate"), type="integer", default=0.75,
                help="intrinsic growth rate"));
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser);
  return(args)
}

args <- parseArgs()
r <- args$growthrate
r
