module fortran_api_test
  use, intrinsic :: iso_c_binding, only : c_bool, c_char, c_double, c_int, c_int32_t, c_int64_t, c_size_t, &
    & c_associated
  use testdrive, only : new_unittest, unittest_type, error_type, check, test_failed
  use gauxc_status, only : gauxc_status_type, gauxc_status_message
  use gauxc_types, only : gauxc_type_molecule, gauxc_type_basisset, gauxc_type_molgrid, &
    & gauxc_type_runtime_environment, gauxc_type_load_balancer, gauxc_type_load_balancer_factory, &
    & gauxc_type_molecular_weights, gauxc_type_molecular_weights_factory, gauxc_type_functional, &
    & gauxc_type_integrator
  use gauxc_enums, only : gauxc_executionspace, gauxc_pruningscheme, gauxc_radialquad, &
    & gauxc_atomicgridsizedefault, gauxc_supportedalg, gauxc_xcweightalg
  use gauxc_atom, only : gauxc_atom_type
  use gauxc_shell, only : gauxc_shell_type
  use gauxc_molecule, only : gauxc_molecule_type, gauxc_molecule_new, gauxc_molecule_delete, &
    & gauxc_molecule_new_from_atoms, gauxc_molecule_natoms, gauxc_molecule_equal
  use gauxc_basisset, only : gauxc_basisset_type, gauxc_basisset_new_from_shells, gauxc_basisset_delete
  use gauxc_molgrid, only : gauxc_molgrid_type, gauxc_molgrid_new_default, gauxc_molgrid_delete
  use gauxc_runtime_environment, only : gauxc_runtime_environment_type, gauxc_runtime_environment_new, &
    & gauxc_runtime_environment_delete, gauxc_runtime_environment_comm_rank, gauxc_runtime_environment_comm_size
  use gauxc_load_balancer, only : gauxc_load_balancer_type, gauxc_load_balancer_factory_type, &
    & gauxc_load_balancer_factory_new, gauxc_load_balancer_factory_delete, &
    & gauxc_load_balancer_factory_get_instance, gauxc_load_balancer_delete
  use gauxc_molecular_weights, only : gauxc_molecular_weights_type, gauxc_molecular_weights_factory_type, &
    & gauxc_molecular_weights_settings, gauxc_molecular_weights_factory_new, &
    & gauxc_molecular_weights_factory_delete, gauxc_molecular_weights_factory_get_instance, &
    & gauxc_molecular_weights_delete, gauxc_molecular_weights_modify_weights
  use gauxc_xc_functional, only : gauxc_functional_type, gauxc_functional, gauxc_functional_from_enum, &
    & gauxc_functional_from_string, gauxc_functional_delete
  use gauxc_integrator, only : gauxc_integrator_type, gauxc_integrator_new, gauxc_integrator_delete
  implicit none
  private

  public :: collect_fortran_api_suite

contains

  subroutine collect_fortran_api_suite(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)

    testsuite = [ &
      new_unittest("status-message", test_status_message), &
      new_unittest("enum-correspondence", test_enum_correspondence), &
      new_unittest("molecule-lifecycle", test_molecule_lifecycle), &
      new_unittest("basisset-header", test_basisset_header), &
      new_unittest("runtime-environment", test_runtime_environment), &
      new_unittest("functional-create-delete", test_functional_create_delete), &
      new_unittest("pipeline-header-tags", test_pipeline_header_tags) &
    ]

  end subroutine collect_fortran_api_suite

  subroutine expect_status_ok(error, status, where)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type), intent(in) :: status
    character(len=*), intent(in) :: where

    call check(error, status%code, 0_c_int)
    if (allocated(error)) then
      call test_failed(error, where, trim(gauxc_status_message(status)))
      return
    end if
  end subroutine expect_status_ok

  subroutine fill_water_atoms(atoms)
    type(gauxc_atom_type), intent(out) :: atoms(3)

    atoms(1)%atomic_number = 1_c_int64_t
    atoms(1)%x = 0.0_c_double
    atoms(1)%y = 1.579252144093028_c_double
    atoms(1)%z = 2.174611055780858_c_double

    atoms(2)%atomic_number = 8_c_int64_t
    atoms(2)%x = 0.0_c_double
    atoms(2)%y = 0.0_c_double
    atoms(2)%z = 0.0_c_double

    atoms(3)%atomic_number = 1_c_int64_t
    atoms(3)%x = 0.0_c_double
    atoms(3)%y = 1.579252144093028_c_double
    atoms(3)%z = -2.174611055780858_c_double
  end subroutine fill_water_atoms

  subroutine fill_minimal_shells(shells)
    type(gauxc_shell_type), intent(out) :: shells(3)

    shells = gauxc_shell_type( &
      l = 0_c_int32_t, pure = .false._c_bool, nprim = 0_c_int32_t, &
      exponents = 0.0_c_double, coefficients = 0.0_c_double, origin = 0.0_c_double, &
      shell_tolerance = -1.0_c_double)

    shells(1)%l = 0_c_int32_t
    shells(1)%pure = .false._c_bool
    shells(1)%nprim = 1_c_int32_t
    shells(1)%exponents(1) = 3.42525091_c_double
    shells(1)%coefficients(1) = 1.0_c_double
    shells(1)%origin = [0.0_c_double, 1.579252144093028_c_double, 2.174611055780858_c_double]
    shells(1)%shell_tolerance = -1.0_c_double

    shells(2)%l = 0_c_int32_t
    shells(2)%pure = .false._c_bool
    shells(2)%nprim = 1_c_int32_t
    shells(2)%exponents(1) = 130.709321_c_double
    shells(2)%coefficients(1) = 1.0_c_double
    shells(2)%origin = [0.0_c_double, 0.0_c_double, 0.0_c_double]
    shells(2)%shell_tolerance = -1.0_c_double

    shells(3)%l = 0_c_int32_t
    shells(3)%pure = .false._c_bool
    shells(3)%nprim = 1_c_int32_t
    shells(3)%exponents(1) = 3.42525091_c_double
    shells(3)%coefficients(1) = 1.0_c_double
    shells(3)%origin = [0.0_c_double, 1.579252144093028_c_double, -2.174611055780858_c_double]
    shells(3)%shell_tolerance = -1.0_c_double
  end subroutine fill_minimal_shells

  subroutine test_status_message(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status

    status%code = 0_c_int
    call check(error, gauxc_status_message(status), "")
    if (allocated(error)) return
  end subroutine test_status_message

  subroutine test_enum_correspondence(error)
    type(error_type), allocatable, intent(out) :: error

    call check(error, gauxc_executionspace%host, 0_c_int)
    if (allocated(error)) return
    call check(error, gauxc_executionspace%device, 1_c_int)
    if (allocated(error)) return

    call check(error, gauxc_supportedalg%xc, 0_c_int)
    if (allocated(error)) return
    call check(error, gauxc_supportedalg%den, 1_c_int)
    if (allocated(error)) return
    call check(error, gauxc_supportedalg%snlink, 2_c_int)
    if (allocated(error)) return

    call check(error, gauxc_pruningscheme%unpruned, 0_c_int)
    if (allocated(error)) return
    call check(error, gauxc_pruningscheme%robust, 1_c_int)
    if (allocated(error)) return
    call check(error, gauxc_pruningscheme%treutler, 2_c_int)
    if (allocated(error)) return
  end subroutine test_enum_correspondence

  subroutine test_molecule_lifecycle(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status
    type(gauxc_atom_type) :: atoms(3)
    type(gauxc_molecule_type) :: mol_a, mol_b, mol_empty
    integer(c_size_t) :: natoms
    logical :: is_equal

    call fill_water_atoms(atoms)

    mol_empty = gauxc_molecule_new(status)
    call expect_status_ok(error, status, "gauxc_molecule_new")
    if (allocated(error)) return
    call check(error, mol_empty%hdr%type, gauxc_type_molecule)
    if (allocated(error)) return

    natoms = gauxc_molecule_natoms(status, mol_empty)
    call expect_status_ok(error, status, "gauxc_molecule_natoms(empty)")
    if (allocated(error)) return
    call check(error, int(natoms, kind=c_int64_t), 0_c_int64_t)
    if (allocated(error)) return

    mol_a = gauxc_molecule_new_from_atoms(status, atoms, int(size(atoms), kind=c_size_t))
    call expect_status_ok(error, status, "gauxc_molecule_new_from_atoms(a)")
    if (allocated(error)) return

    mol_b = gauxc_molecule_new_from_atoms(status, atoms, int(size(atoms), kind=c_size_t))
    call expect_status_ok(error, status, "gauxc_molecule_new_from_atoms(b)")
    if (allocated(error)) return

    natoms = gauxc_molecule_natoms(status, mol_a)
    call expect_status_ok(error, status, "gauxc_molecule_natoms(water)")
    if (allocated(error)) return
    call check(error, int(natoms, kind=c_int64_t), 3_c_int64_t)
    if (allocated(error)) return

    is_equal = gauxc_molecule_equal(status, mol_a, mol_b)
    call check(error, is_equal)
    if (allocated(error)) return
    call expect_status_ok(error, status, "gauxc_molecule_equal")
    if (allocated(error)) return

    call gauxc_molecule_delete(status, mol_b)
    call expect_status_ok(error, status, "gauxc_molecule_delete(b)")
    if (allocated(error)) return

    call gauxc_molecule_delete(status, mol_a)
    call expect_status_ok(error, status, "gauxc_molecule_delete(a)")
    if (allocated(error)) return

    call gauxc_molecule_delete(status, mol_empty)
    call expect_status_ok(error, status, "gauxc_molecule_delete(empty)")
  end subroutine test_molecule_lifecycle

  subroutine test_basisset_header(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status
    type(gauxc_shell_type) :: shells(3)
    type(gauxc_basisset_type) :: basis

    call fill_minimal_shells(shells)
    basis = gauxc_basisset_new_from_shells(status, shells, int(size(shells), kind=c_size_t))
    call expect_status_ok(error, status, "gauxc_basisset_new_from_shells")
    if (allocated(error)) return

    call check(error, basis%hdr%type, gauxc_type_basisset)
    if (allocated(error)) return

    call gauxc_basisset_delete(status, basis)
    call expect_status_ok(error, status, "gauxc_basisset_delete")
  end subroutine test_basisset_header

  subroutine test_runtime_environment(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status
    type(gauxc_runtime_environment_type) :: rt
    integer(c_int) :: rank, nproc

    rt = gauxc_runtime_environment_new(status)
    call expect_status_ok(error, status, "gauxc_runtime_environment_new")
    if (allocated(error)) return

    call check(error, rt%hdr%type, gauxc_type_runtime_environment)
    if (allocated(error)) return

    rank = gauxc_runtime_environment_comm_rank(status, rt)
    call expect_status_ok(error, status, "gauxc_runtime_environment_comm_rank")
    if (allocated(error)) return

    nproc = gauxc_runtime_environment_comm_size(status, rt)
    call expect_status_ok(error, status, "gauxc_runtime_environment_comm_size")
    if (allocated(error)) return

    call check(error, nproc >= 1_c_int)
    if (allocated(error)) return
    call check(error, rank >= 0_c_int)
    if (allocated(error)) return
    call check(error, rank < nproc)
    if (allocated(error)) return

    call gauxc_runtime_environment_delete(status, rt)
    call expect_status_ok(error, status, "gauxc_runtime_environment_delete")
  end subroutine test_runtime_environment

  subroutine test_functional_create_delete(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status
    type(gauxc_functional_type) :: func_from_str, func_from_enum
    character(kind=c_char, len=*), parameter :: svwn5 = "SVWN5"

    func_from_str = gauxc_functional_from_string(status, svwn5, .false._c_bool)
    call expect_status_ok(error, status, "gauxc_functional_from_string")
    if (allocated(error)) return
    call check(error, func_from_str%hdr%type, gauxc_type_functional)
    if (allocated(error)) return

    func_from_enum = gauxc_functional_from_enum(status, gauxc_functional%pbe0, .false._c_bool)
    call expect_status_ok(error, status, "gauxc_functional_from_enum")
    if (allocated(error)) return
    call check(error, func_from_enum%hdr%type, gauxc_type_functional)
    if (allocated(error)) return

    call gauxc_functional_delete(status, func_from_enum)
    call expect_status_ok(error, status, "gauxc_functional_delete(enum)")
    if (allocated(error)) return

    call gauxc_functional_delete(status, func_from_str)
    call expect_status_ok(error, status, "gauxc_functional_delete(string)")
  end subroutine test_functional_create_delete

  subroutine test_pipeline_header_tags(error)
    type(error_type), allocatable, intent(out) :: error
    type(gauxc_status_type) :: status
    type(gauxc_atom_type) :: atoms(3)
    type(gauxc_shell_type) :: shells(3)
    type(gauxc_molecule_type) :: mol
    type(gauxc_basisset_type) :: basis
    type(gauxc_molgrid_type) :: mg
    type(gauxc_runtime_environment_type) :: rt
    type(gauxc_load_balancer_factory_type) :: lbf
    type(gauxc_load_balancer_type) :: lb
    type(gauxc_molecular_weights_settings) :: mw_settings
    type(gauxc_molecular_weights_factory_type) :: mwf
    type(gauxc_molecular_weights_type) :: mw
    type(gauxc_functional_type) :: func
    type(gauxc_integrator_type) :: integ
    character(kind=c_char, len=*), parameter :: input_type = "Replicated"
    character(kind=c_char, len=*), parameter :: kernel_default = "Default"
    character(kind=c_char, len=*), parameter :: func_name = "PBE"

    call fill_water_atoms(atoms)
    call fill_minimal_shells(shells)

    mol = gauxc_molecule_new_from_atoms(status, atoms, int(size(atoms), kind=c_size_t))
    call expect_status_ok(error, status, "pipeline: molecule")
    if (allocated(error)) return
    call check(error, mol%hdr%type, gauxc_type_molecule)
    if (allocated(error)) return
    call check(error, c_associated(mol%ptr), "mol pointer not associated")
    if ((allocated(error))) return

    basis = gauxc_basisset_new_from_shells(status, shells, int(size(shells), kind=c_size_t))
    call expect_status_ok(error, status, "pipeline: basis")
    if (allocated(error)) return
    call check(error, basis%hdr%type, gauxc_type_basisset)
    if (allocated(error)) return
    call check(error, c_associated(basis%ptr), "basis pointer not associated")
    if (allocated(error)) return

    rt = gauxc_runtime_environment_new(status)
    call expect_status_ok(error, status, "pipeline: runtime")
    if (allocated(error)) return
    call check(error, rt%hdr%type, gauxc_type_runtime_environment)
    if (allocated(error)) return
    call check(error, c_associated(rt%ptr), "rt pointer not associated")
    if (allocated(error)) return

    mg = gauxc_molgrid_new_default(status, mol, gauxc_pruningscheme%unpruned, 256_c_int64_t, &
      & gauxc_radialquad%mura_knowles, gauxc_atomicgridsizedefault%finegrid)
    call expect_status_ok(error, status, "pipeline: molgrid")
    if (allocated(error)) return
    call check(error, mg%hdr%type, gauxc_type_molgrid)
    if (allocated(error)) return
    call check(error, c_associated(mg%ptr), "mg pointer not associated")
    if (allocated(error)) return

    lbf = gauxc_load_balancer_factory_new(status, gauxc_executionspace%host)
    call expect_status_ok(error, status, "pipeline: load_balancer_factory")
    if (allocated(error)) return
    call check(error, lbf%hdr%type, gauxc_type_load_balancer_factory)
    if (allocated(error)) return
    call check(error, c_associated(lbf%ptr), "lbf pointer not associated")
    if (allocated(error)) return

    lb = gauxc_load_balancer_factory_get_instance(status, lbf, rt, mol, mg, basis)
    call expect_status_ok(error, status, "pipeline: load_balancer")
    if (allocated(error)) return
    call check(error, lb%hdr%type, gauxc_type_load_balancer)
    if (allocated(error)) return
    call check(error, c_associated(lb%ptr), "lb pointer not associated")
    if (allocated(error)) return

    mwf = gauxc_molecular_weights_factory_new(status, gauxc_executionspace%host)
    call expect_status_ok(error, status, "pipeline: molecular_weights_factory")
    if (allocated(error)) return
    call check(error, mwf%hdr%type, gauxc_type_molecular_weights_factory)
    if (allocated(error)) return
    call check(error, c_associated(mwf%ptr), "mwf pointer not associated")
    if (allocated(error)) return

    mw = gauxc_molecular_weights_factory_get_instance(status, mwf)
    call expect_status_ok(error, status, "pipeline: molecular_weights")
    if (allocated(error)) return
    call check(error, mw%hdr%type, gauxc_type_molecular_weights)
    if (allocated(error)) return
    call check(error, c_associated(mw%ptr), "mw pointer not associated")
    if (allocated(error)) return

    call gauxc_molecular_weights_modify_weights(status, mw, lb)
    call expect_status_ok(error, status, "pipeline: modify_weights")
    if (allocated(error)) return

    func = gauxc_functional_from_string(status, func_name, .false._c_bool)
    call expect_status_ok(error, status, "pipeline: functional")
    if (allocated(error)) return
    call check(error, func%hdr%type, gauxc_type_functional)
    if (allocated(error)) return
    call check(error, c_associated(func%ptr), "functional pointer not associated")
    if (allocated(error)) return

    integ = gauxc_integrator_new(status, func, lb, gauxc_executionspace%host)
    call expect_status_ok(error, status, "pipeline: integrator")
    if (allocated(error)) return
    call check(error, integ%hdr%type, gauxc_type_integrator)
    if (allocated(error)) return
    call check(error, c_associated(integ%ptr), "integrator pointer not associated")
    if (allocated(error)) return

    call gauxc_integrator_delete(status, integ)
    call expect_status_ok(error, status, "cleanup: integrator")
    if (allocated(error)) return

    call gauxc_functional_delete(status, func)
    call expect_status_ok(error, status, "cleanup: functional")
    if (allocated(error)) return

    call gauxc_molecular_weights_delete(status, mw)
    call expect_status_ok(error, status, "cleanup: molecular_weights")
    if (allocated(error)) return

    call gauxc_molecular_weights_factory_delete(status, mwf)
    call expect_status_ok(error, status, "cleanup: molecular_weights_factory")
    if (allocated(error)) return

    call gauxc_load_balancer_delete(status, lb)
    call expect_status_ok(error, status, "cleanup: load_balancer")
    if (allocated(error)) return

    call gauxc_load_balancer_factory_delete(status, lbf)
    call expect_status_ok(error, status, "cleanup: load_balancer_factory")
    if (allocated(error)) return

    call gauxc_molgrid_delete(status, mg)
    call expect_status_ok(error, status, "cleanup: molgrid")
    if (allocated(error)) return

    call gauxc_runtime_environment_delete(status, rt)
    call expect_status_ok(error, status, "cleanup: runtime")
    if (allocated(error)) return

    call gauxc_basisset_delete(status, basis)
    call expect_status_ok(error, status, "cleanup: basis")
    if (allocated(error)) return

    call gauxc_molecule_delete(status, mol)
    call expect_status_ok(error, status, "cleanup: molecule")
  end subroutine test_pipeline_header_tags

end module fortran_api_test
