      program Geom
      implicit none
      character(len=300) :: path,blendRoute,GeomFile,scriptfile,params
      character(len=900) :: auxtxt  
      CHARACTER(len=300) :: arg
      real :: maxRadio
      CALL getcwd(path)
      write(scriptfile,'(a,a)') trim(path),"/scriptJustPrint.py"
!     scriptfile = "/Users/marshall/Desktop/script.py"
      CALL get_command_argument(1, arg)
      IF (LEN_TRIM(arg) .ne. 0) then
      read(arg,*) auxtxt
      print*,"archivo: ", trim(auxtxt)
      CALL getcwd(path)
      write(GeomFile,'(a,a,a)') trim(path),"/",trim(auxtxt)
      print*, GeomFile
!     GeomFile = "/Users/marshall/Desktop/input.blend"
      else
      stop "Especifique archivo .blend"
      end if
      ! ruta al programa
      call chdir ("/Users/marshall/Applications")
      blendRoute = "./blender.app/Contents/MacOS/blender"
      
      params = "-b -P"
      
      write(auxtxt,*) trim(blendRoute)," ",trim(GeomFile) & 
                     ," ",trim(params)," ",trim(scriptfile)
      call system(auxtxt)
      write(6,*)"ending"
      end
