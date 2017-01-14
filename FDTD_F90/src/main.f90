PROGRAM Main
  USE Grid
  USE FDTD 
  USE TFSF
  USE PrintField
  USE CPML 
  IMPLICIT NONE

  INTEGER :: t, dummy_t, FrameSkip;
  CHARACTER(LEN=3) :: CharToInt;
  CALL GRID_Initialize();
  CALL CPML_Initialize();
  CALL TFSF_Initialize(); 
  CALL Print_Initialize(); 

  dummy_t = 100;
  FrameSkip = TotalTime / 100;
  
  DO t = 1, TotalTime
  
    !==========
    ! 1-D Grid
    !==========
    IF(SRC_TYPE .EQ. 0) THEN
      CALL TFSF_Inc();
      Einc0(-2) = Einc0(-2) + Source(t);
      CALL TFSF_UpdateHinc();
      CALL TFSF_UpdateEinc();
    ELSE
      Ez(SizeX/2,SizeY/2) = Ez(SizeX/2,SizeY/2) + Source(t);
    END IF 
       
    !==========
    ! 2-D Grid
    !==========
    CALL UpdateH();
    
    CALL CPML_UpdateHx();
    CALL CPML_UpdateHy();
    
   
    CALL UpdateE();
   
    IF(SRC_TYPE .EQ. 0) THEN
      CALL TFSF_Hupdate();
      CALL TFSF_Eupdate();
    END IF
      
    CALL CPML_UpdateEz(); 
    
    IF(TotalTime <= 100) THEN
      WRITE(CharToInt,"(I3)"),dummy_t
      CALL PrintEz(CharToInt);
      dummy_t = dummy_t + 1;
    ELSE IF(MOD(t,FrameSkip) .EQ. 0) THEN
      WRITE(CharToInt,"(I3)"),dummy_t
      CALL PrintEz(CharToInt);
      dummy_t = dummy_t + 1;
    END IF 
  
  END DO
    
  CALL Print_Limits();

  CALL TFSF_Finalize();
  CALL CPML_Finalize(); 
  CALL GRID_Finalize();
  
END PROGRAM Main



