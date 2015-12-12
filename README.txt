
CIS 560/460 FINAL PROJECT 


-BVH integration and surface area heuristic

	- successfully integrated the BVH tree implementation into the existing
	 Monte Carlo Lighting Estimator
	- BVH uses surface area heuristic
	-slowed down scene loading compared to strictly coordinate based splitting 
	of tree, but sped up rendering process 



- Scene Files 
	- created 3 different scene files, universe.xml, deer.xml, and sky_box.xml. 
	-sky_box.xml contains two dodecahedrons in a space backdrop box
	-deer.xml contains a deer obj with forest backdrop
	-universe.xml contains the solar system with a space backdrop




-Problems
	- Monte Carlo indirect lighting still not completely functional - made 
	some improvements from the last assignment (was not functional in last assignment)
	fixed bugs from last assginments such as: 
		- wi_ret not being set correctly 
		-indirect lighting causing crashes when isx.object_hit = null
	- Mesh's not rendering properly still - regular mesh renders partially, but not in my own scene - sky_box.xml (crashed in last assignment)


