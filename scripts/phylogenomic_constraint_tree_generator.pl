#!/usr/bin/perl

$scriptname=$0; $scriptname =~ s/.+\///g;

if ($#ARGV != 2 || $ARGV[0] eq "-h" || $ARGV[0] eq "-help" || $ARGV[0] eq "--help")
        {print "\nGenerate a constraint tree according to your input\n";
        print "Usage:\t $scriptname <tree_topology> <species_annotation>  <original_ML_unconstraint_tree>\n";
        print "Example file:\n\n";
        print "tree_topology: (POHYD,Bilateria)\n\n";
        print "species_annotation: POHYD=POHYD
Bilateria=((((Aplysia_californica,Crassostrea_gigas),(Capitella_teleta,Helobdella_robusta)),((Daphnia_pulex,Tribolium_castaneum),Ixodes_scapularis)),((Branchiostoma_floridae,Homo_sapiens),(Saccoglossus_kowalevskii,Strongylocentrotus_purpuratus)))\n\n";
        print "original_ML_unconstraint_tree: ((Codosiga_hollandica_17,(((Choanoeca_perplexa_11,Monosiga_brevicollis_mx1),(Salpingoeca_infusionum_12,(Salpingoeca_roanoka_13,Salpingoeca_rosetta))),((((Mylnosiga_fluctuans_19,Salpingoeca_helianthica_18),(Salpingoeca_macrocollata_06,Salpingoeca_punica_03)),(Salpingoeca_qvevrii_09,Salpingoeca_urceolata_04)),(Salpingoeca_dolichothecata_16,(Didymoeca_costata_10,(Diaphanoeca_grandis_RI_01,(Stephanoeca_diplocostata_AUFR,(Acanthoeca_spectabilis_VA_02,'Acanthoeca_sp._10tr')))))))),(((((((((Aiptasia_pallida,Bolocera_tuedia),(Edwardsiella_lineata,Nematostella_vectensis)),((Antipathes_caribbeana,Plumapathes_pennacea),(Montastraea_faveolata,Porites_australiensis))),(Gorgonia_ventalina,Pennatula_rubra)),(((((((Aurelia_aurita,Stomolophus_meleagris),Pelagia_noctiluca),(Atolla_vanhoeffeni,Periphylla_periphylla)),Alatina_alata),((Craspedacusta_sowerbyi,Liriope_tetraphylla),(Hydra_magnipapillata,(Clytia_hemisphaerica,(Hydractinia_polyclina,Nanomia_bijuga))))),Lucernariopsis_campanulata),POHYD)),((((Aplysia_californica,Crassostrea_gigas),(Capitella_teleta,Helobdella_robusta)),((Daphnia_pulex,Tribolium_castaneum),Ixodes_scapularis)),((Branchiostoma_floridae,Homo_sapiens),(Saccoglossus_kowalevskii,Strongylocentrotus_purpuratus)))),Trichoplax_adhaerens),(((((Amphimedon_queenslandica,Petrosia_ficiformis),((Ephydatia_muelleri,Spongilla_lacustris),((Latrunculia_apicalis,Kirkpatrickia_variolosa),Mycale_phyllophila))),(Chondrilla_nucula,(Ircinia_fasciculata,Pleraplysilla_spinifera))),((Aphrocallistes_vastus,((Rossella_fibulata,Sympagella_nux),Euplectella_aspergillum)),Hyalonema_populiferum)),((Clathrina_coriacea,(((Grantia_compressa,(Sycon_ciliatum,Sycon_coactum)),Leuconia_nivea),Leucosolenia_complicata)),((Plakina_jani,Corticium_candelabrum),(Oscarella_carmela,Oscarella_species))))),(((((((Beroe_abyssicola,'Beroe_sp.'),((Mnemiopsis_leidyi,Cestum_veneris),Bolinopsis_infundibulum)),Dryodora_glandiformis),(Hormiphora_californensis,Pleurobrachia_species)),Lampea_pancerina),(Coeloplana_species,Vallicula_multiformis)),Euplokamis_dunlapae)));\n\n";
        print "\n"; exit;
        }


open IN1, "$ARGV[0]";
open IN2, "$ARGV[1]";
open IN3, "$ARGV[2]";
open IN4, "$ARGV[1]";

open OUT, ">constrain.tre";


while (<IN2>) {
	chomp;
	$_ =~ /(\w+)=(\S+)/;
	$my_hash{$1}= $2;

}


while (<IN3>) {
	chomp;
	 while (/\((\w+|'\S+?'),/g) {
         $hash_all{$1} ++;	
}

	while (/,(\w+|'\S+?')\)/g) {
	$hash_all{$1} ++;
}
}

while (<IN4>) {
        chomp;
         while (/\((\w+|'\S+'),/g) {
         $hash_used{$1} ++;
}

        while (/,(\w+|'\S+')\)/g) {
        $hash_used{$1} ++;
}
}

foreach $element (keys %hash_all) {
	if (exists $hash_used{$element}) {

} else {

	push @array, $element;	
}
}



while (<IN1>) {
	chomp;
	for $line (keys %my_hash) {
	$_ =~ s/$line,/$my_hash{$line},/;
       $_ =~ s/$line\)/$my_hash{$line}\)/;
}
       our $core= $_ ;
           
}


foreach $arr (@array) {
	
	our $rest .= ",$arr";

}


print OUT "(${core}${rest});\n"; 

