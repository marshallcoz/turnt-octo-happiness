all: LHS tester


LHS : LHS.for
	gfortran -c LHS.for -O5
	ar rcv libLHS.a *.o
	ranlib libLHS.a

tester : tester.f90 libLHS.a
	gfortran -o tester.out tester.f90 -L. -lLHS

clean :
	rm lib*