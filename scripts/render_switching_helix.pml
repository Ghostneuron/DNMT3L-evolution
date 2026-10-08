load analysis_v2/structure/9MPP.cif, complex
hide everything

select focus, complex and chain N and resi 315-379
show cartoon, focus
color gray80, focus

select switching_helix, complex and chain N and resi 340-354
color orange, switching_helix

select acidic_patch_counterparts, complex and chain N and resi 348+351
show sticks, acidic_patch_counterparts
color magenta, acidic_patch_counterparts

select historical_h358, complex and chain N and resi 357
show sticks, historical_h358
color cyan, historical_h358

bg_color white
set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, 1
orient focus
zoom focus, 5
ray 1800, 1600
png analysis_v2/structure/9MPP_DNMT3L_switching_helix.png, 0, 0, 300
quit
