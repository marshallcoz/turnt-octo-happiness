      program Geom
      character(len=300) :: path,blendRoute,GeomFile,scriptfile,params
      character(len=900) :: auxtxt  
      CHARACTER(len=32) :: arg
      real :: maxRadio
      CALL getcwd(path)
      write(scriptfile,'(a,a)') trim(path),"/script.py"
!     scriptfile = "/Users/marshall/Desktop/script.py"
      ! rutas
      blendRoute = "./blender.app/Contents/MacOS/blender"
      params = "-b -P"
      
!     GeomFile = "/Users/marshall/Desktop/input.blend"
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then
      read(arg,*) auxtxt
      print*,"archivo: ", trim(auxtxt)
      CALL getcwd(path)
      write(GeomFile,'(a,a,a)') trim(path),"/",trim(auxtxt)
      print*, GeomFile
!     GeomFile = "/Users/marshall/Desktop/input.blend"
      else
      stop "(arg1) Especifique archivo .blend en este directorio"
      end if
      
      CALL get_command_argument(2, arg)
      IF (LEN_TRIM(arg) .ne. 0) then
      read(arg,*) maxRadio
      else
      stop "(arg2) Especifique radio Maximo"
      end if
      write(auxtxt,*) trim(blendRoute)," ",trim(GeomFile) & 
                     ," ",trim(params)," ",trim(scriptfile)," --",maxRadio
      call chdir ("/Users/marshall/Applications")
      write(6,'(a,/,a)')"sending:",trim(auxtxt)
      call system(auxtxt)
      write(6,*)"did it"
      end
