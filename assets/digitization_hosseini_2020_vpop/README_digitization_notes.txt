Hosseini 2020 mosunetuzumab Fig. 5 digitization notes

Files:
- hosseini2020_fig5_digitized_IL6_Tcell_green_gray_trajectories_anchor_scaled.csv
  0.25-day samples from panels A-H. Curves: median_green, p5_green, p95_green, max_gray_envelope.
  Green curves are digitized from the visible green lines and anchor-scaled to printed peak labels where the figure prints those labels. max_gray_envelope is the visible upper envelope across gray individual virtual-patient trajectories at each time; it is not guaranteed to be one virtual patient.
- hosseini2020_fig5_digitized_tumor_waterfall_shapes.csv
  Normalized patient-rank samples from panels I-L at 0.002 rank-fraction spacing. The x-axis in Fig. 5 is unlabeled except "Patients"; approximate_virtual_patient_rank_if_N4500 assumes n=4500, consistent with the NHL VPop size stated for Fig. 4 in the same paper.
- hosseini2020_fig5_digitization_overlay_final.png
  Visual QC overlay of extracted/adjusted trajectories and waterfall contours.
- hosseini2020_fig5_digitization_qc_summary.csv and hosseini2020_fig5_waterfall_qc_summary.csv
  Quick maxima/minima and response-fraction checks.

Axis calibration:
- IL6 panels A-D: log10 y-axis; visible ticks 1, 10, 100, 1000 pg/mL and inferred top decade of 10000 pg/mL.
- T-cell panels E-H: linear 0-100% CD8+CD69+ T-cells.
- Waterfall panels I-L: linear y-axis with -100, 0, 100, and 200% grid/tick positions.
- Time panels A-H: 0-42 days.

Limits:
These are raster digitizations of a published figure, not the raw VPop simulations. Nearly vertical post-dose spikes, dashed percentile gaps, line anti-aliasing, and overlapping gray curves impose uncertainty. Use these as approximate calibration/validation targets, not source data.
