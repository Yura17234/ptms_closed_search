# ==============================/ Функция парсит в строке записи UniProt'а название белка /=============================

def smaller_groups(modif):
    # =================
    if 'hydroxyasparagine' in modif or 'Hydroxyasparagine' in modif: return 'Hydroxyasparagine'

    elif 'hydroxyaspartate' in modif or 'Hydroxyaspartate' in modif: return 'Hydroxyaspartate'

    elif 'hydroxyproline' in modif or 'Hydroxyproline' in modif: return 'Hydroxyproline'

    elif 'hydroxylysine' in modif or 'Hydroxylysine' in modif: return 'Hydroxylysine'

    elif 'hydroxyarginine' in modif or 'Hydroxyarginine' in modif: return 'Hydroxyarginine'

    elif 'hydroxyleucine' in modif or 'Hydroxyleucine' in modif: return 'Hydroxyleucine'

    elif 'hydroxyhistidine' in modif or 'Hydroxyhistidine' in modif: return 'Hydroxyhistidine'
    # =================

    elif 'Methionine sulfoxide' in modif or 'Methionine (R)-sulfoxide' in modif: return 'Sulfoxidated methionine'

    # =================
    elif 'acetylmethionine' in modif or 'Acetylmethionine' in modif: return 'Acetylmethionine'

    elif 'acetylalanine' in modif or 'Acetylalanine' in modif: return 'Acetylalanine'

    elif 'acetyllysine' in modif or 'Acetyllysine' in modif: return 'Acetyllysine'

    elif 'acetylproline' in modif or 'Acetylproline' in modif: return 'Acetylproline'

    elif 'acetylcysteine' in modif or 'Acetylcysteine' in modif: return 'Acetylcysteine'

    elif 'acetylvaline' in modif or 'Acetylvaline' in modif: return 'Acetylvaline'

    elif 'acetylthreonine' in modif or 'Acetylthreonine' in modif: return 'Acetylthreonine'

    elif 'acetylserine' in modif or 'Acetylserine' in modif: return 'Acetylserine'

    elif 'acetylglutamate' in modif or 'Acetylglutamate' in modif: return 'Acetylglutamate'

    elif 'acetylglycine' in modif or 'Acetylglycine' in modif: return 'Acetylglycine'

    elif 'acetylaspartate' in modif or 'Acetylaspartate' in modif: return 'Acetylaspartate'
    # =================
    elif 'Phosphoserine' in modif or 'phosphoserine' in modif: return 'Phosphoserine'

    elif 'Phosphotyrosine' in modif or 'phosphotyrosine' in modif: return 'Phosphotyrosine'

    elif 'Phosphothreonine' in modif or 'phosphothreonine' in modif: return 'Phosphothreonine'

    elif 'Phosphohistidine' in modif or 'phosphohistidine' in modif: return 'Phosphohistidine'
    # =================

    elif 'Sulfotyrosine' in modif: return 'Sulfotyrosine'

    elif 'nitrosocysteine' in modif: return 'Nitrosocysteine'

    elif 'cysteinyl cysteine' in modif: return 'Cysteinyl cysteine'

    elif 'glutathionyl cysteine' in modif: return 'Glutathionylcysteine'

    # =================
    elif 'succinyl)cysteine' in modif: return 'Succinylcysteine'

    elif 'succinyllysine' in modif: return 'Succinyllysine'
    # =================

    elif 'AMP-threonine' in modif: return 'AMPthreonine'

    elif 'pyridoxal phosphate)lysine' in modif: return 'Pyridoxal phosphate lysine'

    # # =================

    elif 'carboxyethyl lysine' in modif: return 'Carboxyethyl lysine'

    elif 'bromotyrosine' in modif: return 'Bromotyrosine'

    elif 'malonyllysine' in modif: return 'Malonyllysine'

    elif 'lactoyllysine' in modif: return 'Lactoyllysine'

    elif 'glutaryllysine' in modif: return 'Glutaryllysine'

    elif 'biotinyllysine' in modif: return 'Biotinyllysine'

    elif 'carboxyglutamate' in modif: return 'Carboxyglutamate'

    elif 'carboxylysine' in modif: return 'Carboxylysine'

    elif 'oxoalanine' in modif: return 'Oxoalanine'

    # =================

    elif 'butyryllysine' in modif or 'Butyryllysine' in modif: return 'Butyryllysine'

    # =================
    elif 'Deamidated glutamine' in modif or 'deamidated glutamine' in modif: return 'Deamidated glutamine'

    elif 'Deamidated asparagine' in modif or 'deamidated asparagine' in modif: return 'Deamidated asparagine'
    # =================
    elif 'glutamyl polyglutamate' in modif: return 'Glutamyl polyglutamate'

    elif 'glutamyl dopamine' in modif: return 'Glutamyl dopamine'

    elif 'glutamyl serotonin' in modif: return 'Glutamyl serotonin'

    elif 'glutamyl glycerylphosphorylethanolamine' in modif: return 'Glutamyl glycerylphosphorylethanolamine'

    elif 'glutamyl histamine' in modif: return 'Glutamyl histamine'

    elif 'glutamyl glycine' in modif: return 'Glutamyl glycine'
    # =================
    elif 'trimethyllysine' in modif or 'Trimethyllysine' in modif: return 'Trimethyllysine'

    elif 'trimethylglycine' in modif or 'Trimethylglycine' in modif: return 'Trimethylglycine'

    elif 'trimethylserine' in modif or 'Trimethylserine' in modif: return 'Trimethylserine'

    elif 'trimethylalanine' in modif or 'Trimethylalanine' in modif: return 'Trimethylalanine'
    # =================
    elif 'dimethylarginine' in modif or 'Dimethylarginine' in modif or 'Dimethylated arginine' in modif or 'dimethylated arginine' in modif:
        return 'Dimethylarginine'

    elif 'dimethylasparagine' in modif or 'Dimethylasparagine' in modif: return 'Dimethylasparagine'

    elif 'dimethyllysine' in modif or 'Dimethyllysine' in modif: return 'Dimethyllysine'

    elif 'dimethylserine' in modif or 'Dimethylserine' in modif: return 'Dimethylserine'

    elif 'dimethylproline' in modif or 'Dimethylproline' in modif: return 'Dimethylproline'

    elif 'dimethylglycine' in modif or 'Dimethylglycine' in modif: return 'Dimethylglycine'
    # =================
    elif 'methylarginine' in modif or 'Methylarginine' in modif or 'methylated arginine' in modif: return 'Methylarginine'

    elif 'methyllysine' in modif or 'Methyllysine' in modif or 'methylated lysine' in modif: return 'Methyllysine'

    elif 'methylglutamine' in modif or 'Methylglutamine' in modif: return 'Methylglutamine'

    elif 'methylhistidine' in modif or 'Methylhistidine' in modif: return 'Methylhistidine'

    elif 'methylcysteine' in modif or 'Methylcysteine' in modif: return 'Methylcysteine'

    elif 'methylserine' in modif or 'Methylserine' in modif: return 'Methylserine'

    elif 'methylglycine' in modif or 'Methylglycine' in modif: return 'Methylglycine'

    # =================
    elif 'ADP-ribosyl glutamic' in modif: return 'ADP-ribosylglutamate'

    elif 'ADP-ribosylglycine' in modif: return 'ADP-ribosylglycine'

    elif 'ADP-ribosylthreonine' in modif: return 'ADP-ribosylthreonine'

    elif 'ADP-ribosylserine' in modif: return 'ADP-ribosylserine'

    elif 'ADP-ribosylcysteine' in modif: return 'ADP-ribosylcysteine'

    elif 'ADP-ribosyldiphthamide' in modif: return 'ADP-ribosyldiphthamide'

    elif 'ADP-ribosylarginine' in modif: return 'ADP-ribosylarginine'

    elif 'ADP-ribosylasparagine' in modif: return 'ADP-ribosylasparagine'

    elif 'ADP-ribosyl)lysine' in modif or 'ADP-ribosyllysine' in modif: return 'ADP-ribosyllysine'
    # =================

    elif 'Allysine' in modif: return 'Allysine'

    elif 'glutamyl polyglutamate' in modif: return 'Glutamyl polyglutamate'

    elif 'nitrotyrosine' in modif: return 'Nitrotyrosine'

    elif 'crotonyllysine' in modif: return 'Crotonyllysine'

    elif 'propionyllysine' in modif: return 'Propionyllysine'

    elif 'hydroxybutyryl' in modif or 'hydroxyisobutyryl' in modif: return 'Hydroxyisobutyryl lysine'
    # ====================================
    elif 'Pyrrolidone carboxylic' in modif: return 'Pyrrolidone carboxylic acid'

    elif 'Cysteine sulfinic acid' in modif or 'Cysteine sulfonic acid' in modif or 'Cysteine sulfenic acid' in modif: return 'Cysteine sulfonic acid'

    elif 'Citrulline' in modif: return 'Citrullinearginine'

    elif 'Cysteine methyl ester' in modif: return 'Cysteine methyl ester'

    elif 'Cysteine persulfide' in modif: return 'Cysteine persulfide'

    elif 'Diphthamide' in modif: return 'Diphthamide'

    elif 'Hypusine' in modif: return 'Hypusinelysine'

    elif 'topaquinone' in modif: return 'Topaquinone tyrosine'

    elif 'Arginine amide' in modif: return 'Arginine amide'

    elif 'Glycine amide' in modif: return 'Glycine amide'

    elif 'lipoyllysine' in modif: return 'Lipoyllysine'

    else:
        return None