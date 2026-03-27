# PROMPT_LOG.md

Below please log the exact prompts I give you, numbered, with timestamps.

1. 2026-03-26 8:23PM 
Determine the latest version of the JWST pipeline:
https://github.com/spacetelescope/jwst

Based on the version, report back the micromamba commands we'll need to install that environment. For example, for version 1.19.1:

micromamba create -n jwst_1.19.1 python=3.13
micromamba activate jwst_1.19.1
pip install jwst==1.19.1

Check for micromamba environments I already have installed and tell me which is most recent and whether I already have that version of the JWST pipeline installed.

2. 2026-03-26 8:26PM 
Thanks! Were there any gotchas we should add to my instructions for next time? You didn't find an environment jwst_1.19.1? I thought I'd created that...

3. 2026-03-26 8:30PM 
Yes, please create a new file ~/NIRSpec/wavext/nirspec_wavext_work/notes/PIPELINE_INSTALLATION.md with good instructions for next time. Then go ahead and create that environment.

4. 2026-03-26 8:36PM 
Excellent thank you! Please also keep a complete prompt log including time stamps like this:

~/NIRSpec/wavext/nirspec_wavext_work/notes/
PROMPT_LOG.md
...

5. 2026-03-26 8:40PM 
Thanks now read this file:
https://github.com/spacetelescope/jwst/blob/main/CONTRIBUTING.md
and create a distilled version with the key points we'll need in:
~/NIRSpec/wavext/nirspec_wavext_work/notes/CONTRIBUTING.md

I am planning to work on NIRSpec wavelength extension updates. Should we make the fork now?

Continue adding to the PROMPT_LOG.md
Add instructions to the beginning of PROMPT_LOG.md explaining what to do, as you can see.

6. 2026-03-26 8:41PM 
Thanks! But for PROMPT_LOG.md I just want you to record the exact prompts I give you, numbered, and the timestamps; nothing else please.

7. 2026-03-26 8:46PM 
Great thank you! Now go ahead and fork. Let's call it feature/nirspec_wavelength_extension
Provide details in a new file FORK.md alongside the other notes/
Where will this fork live on my machine? And/or will I be pushing to my Github?

8. 2026-03-26 9:09PM 
Thanks! I'd like to do a lot of work locally before we even think about pushing to the official repository. That being said, pushing to my Github will be useful to track changes as we go. I'm logged into Github in the browser, and I believe I also have the keys in my local environment.

9. 2026-03-26 9:12PM 
Thanks I went ahead and created the fork in my browser:
https://github.com/dancoe/jwst_nirspec_wavext

Please note this in FORK.md
And note I'm using dancoe@gmail.com, *not* my other account OutOfThisWorld

10. 2026-03-26 9:14PM 
Could you try please? I don't seem to have access rights set up in my terminal for some reason. Maybe it's the shell...

11. 2026-03-26 10:00PM 
Excellent thank you! Please generate 2 more files in notes/
* CHANGE_LOG.md – note at the top that this file will contain detailed logs of all changes made – newest on top – including timestamps
* INSTRUCTIONS.md – I'll call this with every prompt, including in future sessions. So for now, just have it say to use the PROMPT_LOG and CHANGE_LOG

12. 2026-03-26 10:48PM 
Thank you! I've added a bit to notes/
Please generate a file notes/NOTES.md that lists and briefly describes all files in that directory

13. 2026-03-26 10:52PM 
Please begin your first investigation: How will we extract the full spectra from the pipeline for various modes: FS, MOS, IFU, (BOTS?). What changes will be required to code and reference files. And then what wavelength ranges should we expand to? Please consult documents described in @[nirspec_wavext_work/notes/NOTES.md], including online resources. Follow @[nirspec_wavext_work/notes/INSTRUCTIONS.md]. Just plan for now without editing code.

14. 2026-03-26 11:00PM 
Thank you! Please save this as notes/IMPLEMENTATION_PLAN.md

In @[nirspec_wavext_work/notes/WAVELENGTH_RANGES.md], please note the nominal wavelength range in one column, followed by the proposed extension in another (e.g., 0.6 – 3.3) instead of separate for Min / Max

15. 2026-03-26 11:16PM 
Thanks! Please update @[nirspec_wavext_work/notes/IMPLEMENTATION_PLAN.md] to say we'll be extending reference files as needed, *not* extrapolating using pipeline code. I've provided links to my CRDS cache, the server, and pipeline notebooks in @[nirspec_wavext_work/notes/PIPELINE_RUNNING.md].

Then yes, please go ahead and update code and generate new local reference files as warranted, keeping detailed notes in @[nirspec_wavext_work/notes/CHANGE_LOG.md]. 

Outside of the pipeline updates, please also plan to download and analyze the JWST 1492 FS data so we can begin calibrating.

16. 2026-03-26 11:25PM 
Thanks! Please generate a new file notes/CALIBRATION.md
Refer to Parlanti et al. (arXiv:2512.14844)

Note they account for 2nd and 3rd order spectra (doesn't work without that)
f(λ) = intrinsic spectrum
S(λ) = observed spectrum
S(λ) ~ k(λ) f(λ) + a(λ) f(λ/2) + b(λ) f(λ/3)
k(λ) = 1st-order contribution: down to ~30% at reddest extended wavelengths
a(λ) = 2nd-order contribution: < ~5%
b(λ) = 3rd-order contribution: < ~1%
Solving for 3 unknown functions (k, a, b) requires 3 sources with both:
f(λ) from standard extraction
S(λ) from extended extraction

Please make detailed notes of the PID 1492 data in a new file notes/DATA.md
What are all the observations, including gratings, exposure times, etc. and how will they be used?

Again, let's focus first on Fixed Slit data. You mentioned enabling IFU wavext. Did you also enable the other modes? We'll prioritize FS. Then MOS, IFU, BOTS.

17. 2026-03-26 11:27PM 
Yes please first move the analysis scripts to a subdirectory analysis/
Then proceed with the data processing and analysis!

18. 2026-03-26 11:32PM 
Please update @[nirspec_wavext_work/notes/DATA.md] by starting with a high-level summary of the high-level (calib 3) products and then a summary of the level 2 files. The target is always the same, so no need to repeat. Also move Usage Strategy to the top before the complete file list.

After that, then yes please proceed as far as you can get this time with data downloading, reduction, and analysis

19. 2026-03-26 11:36PM 
Yes, please make sure you're in the correct environment jwst_1.20.2 that we set up before (see PIPELINE_INSTALLATION.md) then proceed

Make sure to note that in @[nirspec_wavext_work/notes/INSTRUCTIONS.md]

20. 2026-03-26 11:38PM 
Please debug that error and proceed

21. 2026-03-26 11:47PM 
Are you sure it's working now? It looks like you're still getting this error
If we pick this up in another session tomorrow, where should I point you for instructions? @[nirspec_wavext_work/notes/IMPLEMENTATION_PLAN.md]?

22. 2026-03-26 11:49PM 
Thanks and please also make note in @[nirspec_wavext_work/notes/IMPLEMENTATION_PLAN.md] of these downgrade tests and anything else you tried that hasn't worked yet

23. 2026-03-26 11:51PM
Please attempt to unblock us from the errors described in @[nirspec_wavext_work/notes/IMPLEMENTATION_PLAN.md]. Follow @[nirspec_wavext_work/notes/INSTRUCTIONS.md]

24. 2026-03-27 12:05AM
Please reorganize my directories so that I can have two 

jwst_nirspec_wavext –>
https://github.com/dancoe/jwst_nirspec_wavext

notes, analysis –>
https://github.com/dancoe/nirspec_wavext_work
along with@[pipeline/wavelengthrange_extended.asdf] (moved to nirspec_wavext_work/reference_files/)
which can go in a subdirectory reference_files/

Follow @[nirspec_wavext_work/notes/INSTRUCTIONS.md]

25. 2026-03-27 12:03AM
Thank you! Please commit nirspec_wavext_work to my Github dancoe account if you haven't already. And note in @[nirspec_wavext_work/notes/INSTRUCTIONS.md] you should update commit every time we make updates
