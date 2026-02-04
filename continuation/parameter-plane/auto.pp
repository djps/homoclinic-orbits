#==============
# AUTO File PP
#==============

# Load the files pp.f and c.pp.1 into the AUTO command interpreter.
pp = load(equation="pp")
pp = load("pp",constants="pp.1")

# Run and store the result in the Python variable first
first = run("pp")
save(first,"mag.1")

# Set the new start label to the first LP label in b.first and s.first
#second = run(first("UZ1"),e='mag',constants="mag.2")
#save(second,"mag.2")

#mag = load("mag",constants="mag.2")
#second = run(first("UZ1"),e='mag',c='mag.2')
#save(second,"mag.2")

# Set the new start label 
#third = run(first("UZ1"),e='mag',c='mag.3')
#save(third,"mag.3")

#fourth = run(first("UZ3"),e='mag',c='mag.2')
#save(fourth,"mag.4")
