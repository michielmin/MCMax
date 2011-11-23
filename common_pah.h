c
c     Temperature grid
c
      common/temp/tempgrid,dtempgrid,cooltime,qcool
      doubleprecision tempgrid(QUANTUM_NTEMP)
      doubleprecision dtempgrid(QUANTUM_NTEMP)
      doubleprecision cooltime(QUANTUM_NTEMP)
      doubleprecision qcool(QUANTUM_NTEMP)
      common/itemp/ntempgrid
      integer ntempgrid
c
c     Input field
c
      common/inprad/radfield
      doubleprecision radfield(FRSIZE_FREQ)
c
c     Emissivity
c
      common/emis/emissivity,emistot
      doubleprecision emissivity(FRSIZE_FREQ),emistot
c
c     Heat-content grid
c
      common/heat/heatcontent
      doubleprecision heatcontent(QUANTUM_NTEMP)
c
c     Peak temperatures
c
      common/peak/peak_temp,peak_eps
      doubleprecision peak_temp(QUANTUM_NTEMP,FRSIZE_FREQ)
      doubleprecision peak_eps(QUANTUM_NTEMP,FRSIZE_FREQ)
      common/ipeak/peak_itemp
      integer peak_itemp(QUANTUM_NTEMP,FRSIZE_FREQ)
c
c     PAH Destruction rate
c     
      common/pahdestr/pahdes_temp,pahdes_time
      doubleprecision pahdes_temp(DUST_SPECIES_MAX)
      doubleprecision pahdes_time(DUST_SPECIES_MAX)
c
c     Matrix for the temperature distribution function
c
      common/tempdist/mat,tdistrib
      doubleprecision mat(QUANTUM_NTEMP,QUANTUM_NTEMP)
      doubleprecision tdistrib(QUANTUM_NTEMP)
c
c     Information for the user
c
      common/quantinfo/peak_temp_store,cooltime_store
c     %                 pahdes_time_store
      doubleprecision peak_temp_store(QUANTUM_NTEMP,
     %                     FRSIZE_FREQ,DUST_SPECIES_MAX)
      doubleprecision cooltime_store(QUANTUM_NTEMP,DUST_SPECIES_MAX)
c      doubleprecision pahdes_time_store(DUST_SPECIES_MAX,
c     %                            FRSIZE_Y_SMALL,FRSIZE_X)
