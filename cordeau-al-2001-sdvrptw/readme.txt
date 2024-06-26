
=======================================================================================
                    SDVRPTW Instances of Cordeau et al. (1997)

Solution format, as specified at  http://neumann.hec.ca/chairedistributique/data/:
========================================================================================

The format of data and solution files in all directories (except "darp", "gmstp", "gqap") is as follows:

A) DATA FILES

The first line contains the following information:

	type m n t

where

	type = 0 (VRP)
               1 (PVRP)
               2 (MDVRP)
               3 (SDVRP)
               4 (VRPTW)
               5 (PVRPTW)
               6 (MDVRPTW)
               7 (SDVRPTW)

	m = number of vehicles

	n = number of customers

	t = number of days (PVRP), depots (MDVRP) or vehicle types (SDVRP)


The next t lines contain, for each day (or depot or vehicle type), the 
following information:

        D Q

where

	D = maximum duration of a route

	Q = maximum load of a vehicle


The next lines contain, for each customer, the following information:

	i x y d q f a list e l

where

	i = customer number

	x = x coordinate

	y = y coordinate

	d = service duration

	q = demand

	f = frequency of visit

	a = number of possible visit combinations

	list = list of all possible visit combinations

        e = beginning of time window (earliest time for start of service),
            if any

        l = end of time window (latest time for start of service), if any
 
               Each visit combination is coded with the decimal equivalent of
               the corresponding binary bit string. For example, in a 5-day
               period, the code 10 which is equivalent to the bit string 01010
               means that a customer is visited on days 2 and 4. (Days are
               numbered from left to right.)

Note : In the case of the MDVRP, the lines go from 1 to n + t and the last t
entries correspond to the t depots. In the case of the VRP, PVRP and SDVRP, the
lines go from 0 to n and the first entry corresponds to the unique depot.



B) SOLUTION FILES

The first line contains the cost of the solution (total duration excluding
service time).

The next lines contain, for each route, the following information:

	l k d q list

where

	l = number of the day (or depot or vehicle type)

	k = number of the vehicle

	d = duration of the route

	q = load of the vehicle

	list = ordered sequence of customers (with start-of-service times, if
               applicable)

