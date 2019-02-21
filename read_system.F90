
subroutine read_system_data(sys,in_file)
    !Description:
    !!@details
    !>Subroutine read_system_data initialises the molsysdat data type sys
    !
    ! Example input file: system.in
    !
    ! IN_SYSTEM file for sssh
    ! enter the system dimensionality
    ! 3
    ! enter the number of atoms
    ! 4
    !Specify the system Geometry (Element, Mass, Position for each atom on a seperate line)
    ! S   31.9720707   -0.855555411012  -0.21877666      0.00263869
    ! S   31.9720707   -0.05287651       0.40679458     -0.00499009
    ! S   31.9720707    0.87923688      -0.17790629      0.02338821
    ! H   1.0007825     0.88292557      -0.32172745     -0.66934001
    ! The number of atoms in fragments a and b respectively are
    ! 4,0
    ! The atom identification numbers in fragment a are
    ! 1,2,3,4
    ! The atom identification numbers in fragment b are
    !
    ! enter the bond lengths, in order, which give the cutoff for a
    ! real bond to exist between atoms i and j (units are bohr)
    ! 5.0
    ! 5.0
    ! 5.0
    ! 5.0
    ! 5.0
    ! 5.0
    !
    !\param[out] sys The molsysdat object to be initialised
    !\param[in] in_file The relative path specifying the IN_SYSTEM file used to initialis the molsysdat object
    !

    use molecule_specs
    implicit none

    type (molsysdat), intent(out) :: sys
    character(len=80), intent(in) :: in_file

    integer :: n,i,j,k,dimen
    integer :: length, ifin
    character(len=80) :: title,comment_line, elements

    !  read comment, which is a header for the SYSTEM file

    write(11,*) '---------------------------------------------------'
    write(11,*) ''
    write(11,*) '   read in molecule specifications from '//trim(in_file)
    write(11,*) ''

    open(unit=7,file=trim(in_file),status='old')

    read(7,80) title
    write(11,81) title
80  format(a79)
81  format(a80)

    !  read in the dimensionality of the space the system is embedded in
    read(7,80) comment_line
    read(7,*) dimen

    write(11,81) comment_line
    write(11,*) dimen

    !  read in the actual number of atoms

    read(7,80) comment_line
    read(7,*)  n

    write(11,81) comment_line
    write(11,*)  n

    !  create a new molsysdat object

    call new(sys,n,dimen,title)

    !Read in the system geometry

    !  read in equilibrium or near-equilibrium geometry
    read(7,80) comment_line
    write(11,80) comment_line
    do i=1,sys%natom
        read(7,*) sys%atom_label(i), sys%mass(i), (sys%EquilibriumGeom(k,i),k=1,sys%dimen)
        write(11,*) sys%atom_label(i), sys%mass(i), (sys%EquilibriumGeom(k,i),k=1,sys%dimen)
    enddo

    ! convert mass in amu to atomic units
    sys%mass = sys%mass*amu2au

    write(11,*) 'atomic masses in atomic units (which are 1822*amu) '
    write(11,*) (sys%mass(i),i=1,sys%natom)

    !  given that we need all the atom-atom bonds,
    !  we just assign the indirect addresses of the bonds

    k=1
    do i=1,sys%natom-1
        do j=i+1,sys%natom
            sys%mb(k)=i
            sys%nb(k)=j
            k=k+1
        enddo
    enddo

    close(unit=7)

    return
end

