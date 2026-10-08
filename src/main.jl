# Revise, if installed, applies edits to src/ files without restarting the REPL
try using Revise catch end
using MsiaGen


# Sites (folders in data/); the public repository includes Serdang only:
#
#   Peninsular Malaysia                 Sabah and Sarawak
#   --------------------------------    -----------------
#   Alor-Setar         Layang-Layang    Bintulu
#   Banting            Lubok-Merbau     Kota-Kinabalu
#   Kluang             Melaka           Kuching
#   Kota-Bharu         Pagoh            Miri
#   Kuala-Pilah        Paloh            Sandakan
#   Kuala-Terengganu   Serdang          Sibu
#   Kuantan            Sitiawan         Tawau
#                      Teluk-Intan
#                      Temerloh

site = "Serdang"     # select a site from the list above
seed = -1           # -1 = random seed; or a fixed number to repeat a run
use_stats = true    # use site statistics to generate weather (true),
                    # or false to generate stats from observed weather

# plot thresholds (change to suit the crop):
dry_day = 0.5       # a day with rain below this (mm) is a dry day
dry_spell = 14      # dry spells of at least this many days are shown
dry_month = 100.0   # a month with rain below this (mm) is a dry month
hot_day = 33.0      # a day with tmax at or above this (°C) is a hot day

folder = normpath(@__DIR__, "..", "data")
thresholds = (; dry_day, dry_spell, dry_month, hot_day)

generate_weather(site; folder=folder, from_stats=use_stats, seed=seed)
plot_weather(site; folder=folder, thresholds...)
check_fit(site; folder=folder)

;
