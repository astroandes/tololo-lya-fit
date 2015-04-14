first_guess.png: analytic_solution.x mc.py 
	python mc.py
	rm model.dat

analytic_solution.x: analytic_solution.c
	cc analytic_solution.c -o analytic_solution.x -lm

clean:
	rm *.x *.dat *.png
