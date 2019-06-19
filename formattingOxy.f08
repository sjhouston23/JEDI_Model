module formattingOxy
implicit none

!****************************** Header Variables *******************************
character(len=*),parameter :: &
H01="(' H^+ production rate [cm^-3 s^-1] from oxygen ion precipitation.',/,&
      ' Input of 1 ion/cm^2/s.')",&
H02="(' H_2^+ production rate [cm^-3 s^-1] from oxygen ion precipitation.',/,&
      ' Input of 1 ion/cm^2/s.')",&
H03="(' H_2 excitation rate [cm^-3 s^-1] from oxygen ion precipitation.',/,&
      ' Input of 1 ion/cm^2/s.',/,&
      ' H_2 is only excited from a single NSIM process, TEX.',/,&
      ' Alt [km] Production Rate')",&
H04="(' Charge state equilibrium fractions.',/,' Energy [keV/u] ONeg',6x,'O',&
      9x,'O+',8x,'O2+',7x,'O3+',7x,'O4+',7x,'O5+',7x,'O6+',7x,'O7+',7x,'O8+')",&
H05="(' E',14x,'SP        Sig       dE     dN        SP1           Ions',/,&
      ' Various data to calculate stopping power for input of 1 ion/cm^2/s.',/,&
      ' Conversion from [eV*cm^2] to [(MeV/mg)*cm^2] is 3.345e-15.',/,&
      ' Ion Energy     dE/dN     Sigma       dE     dN        Sigma*dE     ',&
      ' Counts',/,' [keV/u]        [eV*cm^2] [cm^2]      [eV]   [cm^-2]  ',&
      ' [eV*cm^2]    per bin')",&
H06="(' Alt [km]   ONeg',7x,'O',10x,'O+',9x,'O2+',8x,'O3+',8x,'O4+',8x,'O5+',&
      8x,'O6+',8x,'O7+',8x,'O8+')",&
H07="(32x,'SS           DS         SPEX         DPEX          Sum')",&
H08="(' Altitude integrated photon production [photons cm^-2 s^-1] as a',&
      ' function of charge state.')",&
H09="(' ΔAlt [km]  ONeg',7x,'O',10x,'O+',9x,'O2+',8x,'O3+',8x,'O4+',8x,'O5+',&
      8x,'O6+',8x,'O7+',8x,'O8+')",&
H10="(' Photon production [photons cm^-3 s^-1] as a function of altitude',&
      ' and charge state.')",&
H11="(' Initial input of 1 ion cm^-2 s^-1')",&
!******************************* Notes Variables *******************************
N01="(' Note: The photon production below is from direct excitation of a',&
      ' charge state. If an electron gets excited and then',/,6x,&
      ' relaxes back down to a more energy favorable state, it will emit a',&
      ' photon. Only O6+ and O7+ are capable of',/,6x,&
      ' producing X-Rays. The charge states below will produce in the UV',&
      ' spectrum. O8+ is meaningless here - it is only',/,6x,&
      ' output as a check that the model is outputting correctly. ONeg should',&
      ' also always be 0.00E+00.')",&
N02="(' Note: These charge states are the resultant charge state from a',&
      ' collision. E.g. the production under O7+ is the number',/,6x,&
      ' of photons produced from O8+ gaining an electron that cascades',&
      ' resulting in O7+. Therefore, only O6+',/,6x,&
      ' and O7+ are capable of producing X-Rays. The charge states below will',&
      ' produce in the UV spectrum. O8+ is,',/,6x,&
      ' meaningless here - it is only output as a check that the model is',&
      ' outputting correctly. ONeg should also always,',/,6x,&
      ' be 0.00E+00.')",&
!************************** Data Formatting Variables **************************
!* Notes:
!*   1x to create a space before the data
!*   F7.2 is used for atmosphere ==> "3000.00" or " -88.00"
!*   F8.2 is used for energy ==> "25000.00" or "  150.00"
!*   ES9.2 is used for large/small data values ==> "-1.34E+04" or " 2.31E-11"
!*   Electron two-stream formatting - 912 (F2Str)
!*******************************************************************************
F01="(1x,F7.2,36(2x,ES8.2),2x,ES8.2)",&
F02="(1x,F7.2,2x,ES8.2)",&
F03="(1x,F8.2,5x,10(2x,ES8.2))",&
F04="(1x,F8.2,6x,2(ES9.2,1x),F8.2,2(1x,ES9.2),1x,I11)",&
F05="(1x,F8.2,10(2x,ES9.2))",&
F06="(1x,A8,6(I12,1x))",&
F07="(1x,A8,6(F12.2,1x))",&
F2Str="(1x,1P10E11.3)"  !Electron 2-stream formatting (912)
!****************************** Processes Header *******************************
character(len=10) HProc(36)
character(len=8) Coll(8)
data HProc/&
'DC+DS    ','DC+DPEX  ','DCAI+DS  ','DCAI+DPEX','TI+DS    ','TI+DPEX  ',&
'SC+DS    ','SC+DPEX  ','DC+SS    ','DC+SPEX  ','DCAI+SS  ','DCAI+SPEX',&
'TI+SS    ','TI+SPEX  ','TEX+DS   ','TEX+DPEX ','SI+DS    ','SI+DPEX  ',&
'SC+SS    ','SC+SPEX  ','DI+DS    ','DI+DPEX  ','DCAI     ','DC       ',&
'TI       ','TEX+SS   ','TEX+SPEX ','SI+SS    ','SI+SPEX  ','SC       ',&
'DI+SS    ','DI+SPEX  ','TEX      ','SI       ','DI       ','NEG      '/
data Coll/'NEG     ','SI      ','DI      ','TI      ','SC      ','DC      ',&
          'DCAI    ','TEX     '/
end module formattingOxy
