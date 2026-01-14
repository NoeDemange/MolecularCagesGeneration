# Compiler and flags
CC = gcc
CFLAGS = -O2 -DWEIGHTED

# Executable names
GEN_ICT_EXEC = Gen_ICT
# Source files
GEN_ICT_SRC = Gen_ICT.c

# Action for building both executables
all: clean $(GEN_ICT_EXEC) run_gen_ict

# Compile Gen_ICT
$(GEN_ICT_EXEC): $(GEN_ICT_SRC)
	$(CC) $(CFLAGS) -o $(GEN_ICT_EXEC) $(GEN_ICT_SRC)

# Run Gen_ICT (with input_file as argument)
run_gen_ict: $(GEN_ICT_EXEC)
	./$(GEN_ICT_EXEC) tests/realYARZUN03_YOHPIU_YOHPOA.txt

# Clean up generated files
clean:
	rm -f $(GEN_ICT_EXEC)
