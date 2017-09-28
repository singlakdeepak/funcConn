Node: fast (fsl)
================

 Hierarchy : registration.fast
 Exec ID : fast

Original Inputs
---------------

* args : <undefined>
* bias_iters : <undefined>
* bias_lowpass : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* hyper : <undefined>
* ignore_exception : False
* img_type : <undefined>
* in_files : ['/home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/stripper/MNI152_T1_2mm_brain_brain.nii.gz']
* init_seg_smooth : <undefined>
* init_transform : <undefined>
* iters_afterbias : <undefined>
* manual_seg : <undefined>
* mixel_smooth : <undefined>
* no_bias : <undefined>
* no_pve : <undefined>
* number_classes : <undefined>
* other_priors : <undefined>
* out_basename : <undefined>
* output_biascorrected : <undefined>
* output_biasfield : <undefined>
* output_type : NIFTI_GZ
* probability_maps : <undefined>
* segment_iters : <undefined>
* segments : <undefined>
* terminal_output : stream
* use_priors : <undefined>
* verbose : <undefined>

Execution Inputs
----------------

* args : <undefined>
* bias_iters : <undefined>
* bias_lowpass : <undefined>
* environ : {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
* hyper : <undefined>
* ignore_exception : False
* img_type : <undefined>
* in_files : ['/home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/fast/MNI152_T1_2mm_brain_brain.nii.gz']
* init_seg_smooth : <undefined>
* init_transform : <undefined>
* iters_afterbias : <undefined>
* manual_seg : <undefined>
* mixel_smooth : <undefined>
* no_bias : <undefined>
* no_pve : <undefined>
* number_classes : <undefined>
* other_priors : <undefined>
* out_basename : <undefined>
* output_biascorrected : <undefined>
* output_biasfield : <undefined>
* output_type : NIFTI_GZ
* probability_maps : <undefined>
* segment_iters : <undefined>
* segments : <undefined>
* terminal_output : stream
* use_priors : <undefined>
* verbose : <undefined>

Execution Outputs
-----------------

* bias_field : <undefined>
* mixeltype : <undefined>
* partial_volume_files : ['/home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/fast/MNI152_T1_2mm_brain_brain_pve_0.nii.gz', '/home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/fast/MNI152_T1_2mm_brain_brain_pve_1.nii.gz', '/home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/fast/MNI152_T1_2mm_brain_brain_pve_2.nii.gz']
* partial_volume_map : <undefined>
* probability_maps : <undefined>
* restored_image : <undefined>
* tissue_class_files : <undefined>
* tissue_class_map : <undefined>

Runtime info
------------

* command : fast -S 1 /home/deepak/Desktop/funcConn/funcConnGUI/Scripts/tmp/registration/fast/MNI152_T1_2mm_brain_brain.nii.gz
* duration : 25.80451
* hostname : deepak-OptiPlex-9030-AIO

Terminal output
~~~~~~~~~~~~~~~



Environment
~~~~~~~~~~~

* CLICOLOR : 1
* CLUTTER_IM_MODULE : xim
* COMPIZ_CONFIG_PROFILE : ubuntu
* DBUS_SESSION_BUS_ADDRESS : unix:abstract=/tmp/dbus-5sGUiNPpjZ
* DEFAULTS_PATH : /usr/share/gconf/ubuntu.default.path
* DESKTOP_SESSION : ubuntu
* DISPLAY : :1
* FSLBROWSER : /etc/alternatives/x-www-browser
* FSLDIR : /usr/share/fsl/5.0
* FSLLOCKDIR : 
* FSLMACHINELIST : 
* FSLMULTIFILEQUIT : TRUE
* FSLOUTPUTTYPE : NIFTI_GZ
* FSLREMOTECALL : 
* FSLTCLSH : /usr/bin/tclsh
* FSLWISH : /usr/bin/wish
* GDMSESSION : ubuntu
* GDM_LANG : en_US
* GIT_PAGER : cat
* GNOME_DESKTOP_SESSION_ID : this-is-deprecated
* GNOME_KEYRING_CONTROL : 
* GNOME_KEYRING_PID : 
* GPG_AGENT_INFO : /home/deepak/.gnupg/S.gpg-agent:0:1
* GTK2_MODULES : overlay-scrollbar
* GTK_IM_MODULE : ibus
* GTK_MODULES : gail:atk-bridge:unity-gtk-module
* HOME : /home/deepak
* IM_CONFIG_PHASE : 1
* INSTANCE : Unity
* JOB : gnome-session
* JPY_PARENT_PID : 17516
* LANG : en_IN
* LANGUAGE : en_IN:en
* LD_LIBRARY_PATH : /usr/lib/fsl/5.0
* LESSCLOSE : /usr/bin/lesspipe %s %s
* LESSOPEN : | /usr/bin/lesspipe %s
* LOGNAME : deepak
* LS_COLORS : rs=0:di=01;34:ln=01;36:mh=00:pi=40;33:so=01;35:do=01;35:bd=40;33;01:cd=40;33;01:or=40;31;01:mi=00:su=37;41:sg=30;43:ca=30;41:tw=30;42:ow=34;42:st=37;44:ex=01;32:*.tar=01;31:*.tgz=01;31:*.arc=01;31:*.arj=01;31:*.taz=01;31:*.lha=01;31:*.lz4=01;31:*.lzh=01;31:*.lzma=01;31:*.tlz=01;31:*.txz=01;31:*.tzo=01;31:*.t7z=01;31:*.zip=01;31:*.z=01;31:*.Z=01;31:*.dz=01;31:*.gz=01;31:*.lrz=01;31:*.lz=01;31:*.lzo=01;31:*.xz=01;31:*.bz2=01;31:*.bz=01;31:*.tbz=01;31:*.tbz2=01;31:*.tz=01;31:*.deb=01;31:*.rpm=01;31:*.jar=01;31:*.war=01;31:*.ear=01;31:*.sar=01;31:*.rar=01;31:*.alz=01;31:*.ace=01;31:*.zoo=01;31:*.cpio=01;31:*.7z=01;31:*.rz=01;31:*.cab=01;31:*.jpg=01;35:*.jpeg=01;35:*.gif=01;35:*.bmp=01;35:*.pbm=01;35:*.pgm=01;35:*.ppm=01;35:*.tga=01;35:*.xbm=01;35:*.xpm=01;35:*.tif=01;35:*.tiff=01;35:*.png=01;35:*.svg=01;35:*.svgz=01;35:*.mng=01;35:*.pcx=01;35:*.mov=01;35:*.mpg=01;35:*.mpeg=01;35:*.m2v=01;35:*.mkv=01;35:*.webm=01;35:*.ogm=01;35:*.mp4=01;35:*.m4v=01;35:*.mp4v=01;35:*.vob=01;35:*.qt=01;35:*.nuv=01;35:*.wmv=01;35:*.asf=01;35:*.rm=01;35:*.rmvb=01;35:*.flc=01;35:*.avi=01;35:*.fli=01;35:*.flv=01;35:*.gl=01;35:*.dl=01;35:*.xcf=01;35:*.xwd=01;35:*.yuv=01;35:*.cgm=01;35:*.emf=01;35:*.ogv=01;35:*.ogx=01;35:*.aac=00;36:*.au=00;36:*.flac=00;36:*.m4a=00;36:*.mid=00;36:*.midi=00;36:*.mka=00;36:*.mp3=00;36:*.mpc=00;36:*.ogg=00;36:*.ra=00;36:*.wav=00;36:*.oga=00;36:*.opus=00;36:*.spx=00;36:*.xspf=00;36:
* MANDATORY_PATH : /usr/share/gconf/ubuntu.mandatory.path
* MPLBACKEND : module://ipykernel.pylab.backend_inline
* OLDPWD : /home/deepak/Desktop/funcConn/funcConnGUI
* PAGER : cat
* PATH : /home/deepak/anaconda3/bin:/usr/share/fsl/5.0/bin:/usr/lib/fsl/5.0:/home/deepak/anaconda3/bin:/usr/share/fsl/5.0/bin:/home/deepak/anaconda3/bin:/home/deepak/bin:/home/deepak/.local/bin:/home/deepak/Desktop/funcConn/:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games:/usr/local/games:/snap/bin
* POSSUMDIR : /usr/share/fsl/5.0
* PWD : /home/deepak/Desktop/funcConn/funcConnGUI/Scripts
* QT4_IM_MODULE : xim
* QT_ACCESSIBILITY : 1
* QT_IM_MODULE : ibus
* QT_LINUX_ACCESSIBILITY_ALWAYS_ON : 1
* QT_QPA_PLATFORMTHEME : appmenu-qt5
* SESSIONTYPE : gnome-session
* SHELL : /bin/bash
* SHLVL : 1
* SSH_AUTH_SOCK : /run/user/1000/keyring/ssh
* TERM : xterm-color
* UPSTART_EVENTS : started starting
* UPSTART_INSTANCE : 
* UPSTART_JOB : unity-settings-daemon
* UPSTART_SESSION : unix:abstract=/com/ubuntu/upstart-session/1000/1288
* USER : deepak
* VTE_VERSION : 4205
* WINDOWID : 67153949
* XAUTHORITY : /home/deepak/.Xauthority
* XDG_CONFIG_DIRS : /etc/xdg/xdg-ubuntu:/usr/share/upstart/xdg:/etc/xdg
* XDG_CURRENT_DESKTOP : Unity
* XDG_DATA_DIRS : /usr/share/ubuntu:/usr/share/gnome:/usr/local/share/:/usr/share/:/var/lib/snapd/desktop
* XDG_GREETER_DATA_DIR : /var/lib/lightdm-data/deepak
* XDG_RUNTIME_DIR : /run/user/1000
* XDG_SEAT : seat0
* XDG_SEAT_PATH : /org/freedesktop/DisplayManager/Seat0
* XDG_SESSION_DESKTOP : ubuntu
* XDG_SESSION_ID : c2
* XDG_SESSION_PATH : /org/freedesktop/DisplayManager/Session0
* XDG_SESSION_TYPE : x11
* XDG_VTNR : 7
* XMODIFIERS : @im=ibus
* _ : /home/deepak/anaconda3/bin/jupyter

