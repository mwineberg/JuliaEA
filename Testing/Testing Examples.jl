# Setup one_max population
one_max_parms = OneMaxParameters(32)
bp = Population(24, one_max_parms)
print(bp)
one_max_eval(bp, one_max_parms)
print(bp)

# Setup rastrigin population
rastrigin_parms = RastriginParameters(3);
rGT = Population(24, rastrigin_parms);
rPhT = @time binary2decimal(rGT, rastrigin_parms);
rPhTc = @time binary2decimal_comprehension(rGT, rastrigin_parms);
rPhT2 = decimal2float(rPhT1, rastrigin_parms)