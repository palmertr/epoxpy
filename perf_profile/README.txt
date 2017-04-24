The .prof files are performance profile data collected from epoxpy runs using the CProfile tool as shown below:

'python -m cProfile -o freud_bonding.prof freud_bonding.py'

The profile can then be visualized using tools such as snakeviz as shown below:

'snakeviz freud_bonding.prof'
