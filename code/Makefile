parameters.png: sampler_flatchain.dat plots.py
	python plots.py

sampler_flatchain.dat: analytic_solution.x mc.py
	time python mc.py
	rm model.dat

analytic_solution.x: analytic_solution.c
	cc analytic_solution.c -o analytic_solution.x -lm

clean:
	rm *.x *.dat *.png
