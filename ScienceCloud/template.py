#! /usr/bin/env python

import os
from os.path import abspath, basename
import sys

from gc3libs import Application
from gc3libs.cmdline import SessionBasedScript, \
	existing_file, positive_int
from gc3libs.quantity import GB
	
#
# 1) Run script
# 2) Collect data (manually)

if __name__ == '__main__':
    from template import RScript
    RScript().run()

class RScript(SessionBasedScript):
    """
    This script executes a R-script and produces a .csv file

    Call this python script as:

	    python template.py Rscript.R ~/my-csv-file.csv


    """
    def __init__(self):
        super(RScript, self).__init__(version='1.0')

    def setup_args(self):
	# files
	self.add_param('rscript',type=existing_file,help="R script file to be executed on the remote machine")
	self.add_param('csv_input',type=existing_file,help="A csv file for the script")

    def new_tasks(self, extra):
	# initialize list
        apps_to_run = []
	
	# get full path of the R script file
	r_file = abspath(self.params.rscript)
	csv_file = abspath(self.params.csv_input)
	print r_file
	print csv_file
	# append tasks to the list
	apps_to_run.append(RApp(r_file, csv_file))

	return apps_to_run

class RApp(Application):
    """Run Code"""
    def __init__(self, rscript, csvfile):

	r_script = basename(rscript)
	csv_file = basename(csvfile)

        Application.__init__(
            self,
            arguments=["/usr/bin/Rscript", r_script , csv_file],
            inputs=[rscript,csvfile],
            outputs=["results.csv"],
            output_dir = "code-results",
            stdout="stdout.txt",
            stderr="stderr.txt",
		requested_memory=1*GB)



