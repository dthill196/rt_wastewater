High-throughput files
================

## Files and folder explanations

### Files

#### county_list.txt

This is a text file that is used to iterate through each county. The
`.sub` file initiates a new process for each county.

The `rt_episwer_process.R` script runs for each county.

#### Rt_data_county.rds

This is a file with all the county case and wastewater data that the
process filters for each county and iterates through.

#### rt_episewer.sub

This is the submission file that sends the `rt_episewer_process.R` file
through HTCondor’s high-throughput process.

#### rt_episewer_process.R

This is the R script that performs the analysis. In this example, the
`EpiSewer` Rt model is estimated.

#### rt_episewer_wrapper.sh

This is a wrapper file that selects the process R script for HTCondor to
run.

### Folders

#### output

If you open the `rt_episwer.sub` file in a text editor, you will see
that there is a line directing all output files to be saved to this
`output` folder.

#### logs

If you open the `rt_episwer.sub` file in a text editor, you will see
that there is a line directing all log files to be saved to this `log`
folder.

#### result

If you open the `rt_episwer.sub` file in a text editor, you will see
that there is a line directing all results printed in the console to be
saved to this `result` folder.
