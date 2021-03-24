proc takePicture {myFileNoSuffix} {
    mol new ../results/raw/${myFileNoSuffix}.xyz type xyz first -2 last -1 step 1 filebonds 1 autobonds 1 waitfor all
    mol delrep 0 top
    mol representation CPK 1.000000 0.600000 12.000000 12.000000
    mol color Name
    mol selection {z > -0.9}
    mol material Opaque
    mol addrep top

    mol representation CPK 1.000000 0.600000 12.000000 12.000000
    mol color Name
    mol selection {all}
    mol material Transparent
    mol addrep top

    display projection   Orthographic
    axes location Off
    scale by 1.200000
    scale by 1.200000

    proc vmdrestoremycolors {} {
        set colorcmds {
            {color Display {Background} white}
            {color Display {BackgroundTop} black}
            {color Display {BackgroundBot} blue2}
            {color Display {FPS} white}
            {color Name {C} black}
            {color Name {S} cyan}
            }
            foreach colcmd $colorcmds {
            set val [catch {eval $colcmd}]
          }
    }
    vmdrestoremycolors
    label textsize 1.0

    render snapshot ../results/images/${myFileNoSuffix}_top.bmp

    rotate x by -90
    render snapshot ../results/images/${myFileNoSuffix}_side.bmp
    
    #mol delete top
 }
 
proc analyze {muSi muC} {
    for {set seed 1} {$seed <= 3} {incr seed} {	 
        set file muSi${muSi}_muC${muC}_seed${seed}
        takePicture $file
    }
}

analyze -5.80 -5.55
analyze -5.80 -5.65
#analyze -5.80 -5.75
analyze -5.80 -5.85
analyze -5.80 -5.95
#C:\Users\Will Watkins\Documents\Research Fall 2020 Rappe\aiGCMC\gcmc\plot\muSi-5.80_muC-5.75_seed2_ase.xyz}