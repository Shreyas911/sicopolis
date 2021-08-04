#!/usr/bin/env ./libs/bats/bin/bats

cd ../src

LISDIR="/lis-2.0.30/installation"
NETCDF_FORTRAN_DIR="/usr"
@test "verify that the adjoint code is compiling for GRL" {
	run make -f MakefileTapenade clean; make -f MakefileTapenade driveradjoint HEADER=v5_grl20_ss25ka DOMAIN_SHORT=grl LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
	[ "$status" -eq 0 ]
}

@test "verify that the grdchk code is compiling for GRL" {
        run make -f MakefileTapenade clean; make -f MakefileTapenade drivergrdchk HEADER=v5_grl20_ss25ka DOMAIN_SHORT=grl LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
        [ "$status" -eq 0 ]
}

#@test "verify that the results are correct for GRL" {
#        run make -f MakefileTapenade clean; make -f MakefileTapenade drivergrdchk HEADER=v5_grl20_ss25ka DOMAIN_SHORT=grl LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
#        [ "$status" -eq 0 ]
#
#	run ./drivergrdchk
#	[ "$status" -eq 0 ]
#
#        run make -f MakefileTapenade clean; make -f MakefileTapenade driveradjoint HEADER=v5_grl20_ss25ka DOMAIN_SHORT=grl LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
#        [ "$status" -eq 0 ]
#
#	run ./driveradjoint
#	[ "$status" -eq 0 ]
#
#	cd ../test_ad/
#	run ./compare_ad_fd.sh GradientVals_v5_grl20_ss25ka.dat ../src/AdjointVals_v5_grl20_ss25ka.dat adjoint
#	[ "$status" -eq 0 ]
#}

@test "verify that the adjoint code is compiling for ANT" {
        run make -f MakefileTapenade clean; make -f MakefileTapenade driveradjoint HEADER=v5_ant64_b2_future09_ctrl DOMAIN_SHORT=ant LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
        [ "$status" -eq 0 ]
}

@test "verify that the grdchk code is compiling for ANT" {
        run make -f MakefileTapenade clean; make -f MakefileTapenade drivergrdchk HEADER=v5_ant64_b2_future09_ctrl DOMAIN_SHORT=ant LISDIR=${LISDIR} NETCDF_FORTRAN_DIR=${NETCDF_FORTRAN_DIR} TRAVIS_CI=yes;
        [ "$status" -eq 0 ]
}



