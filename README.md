Roller
======
Proposed Project File Structure:

Roller/
	__init__.py
	Roller.py
	run_roller_example.py
	run_emt.py
	run_goldbetter.py
	
	unittests/
		test_roller.py
		
	data/
		dream4/
			...
			gold_standard.csv
		goldbetter/
			...
			gold_standard.csv
		time_varying/
		emt/
	
	util/
		__init__.py
		Grapher.py
		Evaluator.py
		Permuter.py
		LinearWrapper.py
		ForestWrapper.py
		PLSWrapper.py
	
	output/
	
	
	
Goldbetter/
	__init__.py
	GoldBetter.py

TimeVarying/
	__init__.py
	TimeVarying.py
	
All scripts that run the main process flow will be placed in the base folder in the format "run_...py"
All scripts that output files will print to your current working directory. During development, this will be pointed to output/, which will be ignored by git (.gitignore) during syncing.


