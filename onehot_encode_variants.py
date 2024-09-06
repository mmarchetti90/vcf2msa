#!/usr/bin/env python3

"""
Class to one-hot encode variants given a list of vcf files
"""

### ---------------------------------------- ###

class onehot_variants:
    
    def __init__(self):
        
        pass
    
    ### ------------------------------------ ###
    ### PROCESS VCF FILES                    ###
    ### ------------------------------------ ###
    
    def process_vcfs(self, vcf_files, vcf_ids):
        
        # Store vcf file ids
        self.vcf_ids = vcf_ids
        
        # Init dictionary for variants and their info
        variants_info = {'id' : [],
                         'type' : [],
                         'chr' : [],
                         'pos' : [],
                         'ref' : [],
                         'alt' : [],
                         'ann' : []}
        
        # Init row and col lists (variants and samples, respectively) for creating sparse matrix
        rows, cols = [], []
        
        # Process vcf files
        vcf_header = ['chr', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'info', 'format', 'sample']
        for n_vcf, vcf in enumerate(vcf_files):
            
            if (n_vcf + 1) % 10 == 0 and n_vcf > 0:

                print(f'Processed {n_vcf + 1} / {len(vcf_files)} vcf files', end='\r')
            
            data = pd.read_csv(vcf, sep='\t', comment='#', header=None, names=vcf_header)
            
            # Filter
            data = data.loc[data['filter'].isin(['PASS', 'Pass', 'pass', '.']),]
            
            # Define variants ids
            data['var_id'] = data.chr + '_' + data.pos.astype(str) + '_' + data.ref + '_' + data.alt
            
            # Define variant type
            data['var_type'] = ['substitution' if len(ref) == len(alt) else
                                'deletion' if len(ref) > len(alt) else
                                'insertion' if len(ref) < len(alt) else
                                ''
                                for _,(ref,alt) in data.loc[:, ['ref', 'alt']].iterrows()]
            
            # Get annotation
            data['ann'] = [[i for i in info.split(';') if i.startswith('ANN=')][0] if 'ANN=' in info else
                           ''
                           for info in data['info'].values]
            
            # Find new variants, then add them to variants_info
            new_vars = data.loc[~ data.var_id.isin(variants_info['id']), ]
            variants_info['id'].extend(new_vars.var_id.values)
            variants_info['type'].extend(new_vars.var_type.values)
            variants_info['chr'].extend(new_vars.chr.values)
            variants_info['pos'].extend(new_vars.pos.values)
            variants_info['ref'].extend(new_vars.ref.values)
            variants_info['alt'].extend(new_vars.alt.values)
            variants_info['ann'].extend(new_vars.ann.values)
            
            # Update rows and cols lists
            rows.extend([variants_info['id'].index(var_id) for var_id in data.var_id.values])
            cols.extend([n_vcf for _ in range(data.shape[0])])
        
        # Convert variants_info to pandas data frame
        self.variants_info = pd.DataFrame(variants_info)
        
        # Create sparse matrix
        self.variants_matrix = csc_matrix((np.ones(len(rows)), (rows, cols)), dtype=int)
    
    ### ------------------------------------ ###
    
    def save_data(self, file_prefix='./vcf_data', full_matrix=False):
        
        # Save vcf_ids
        with open(f'{file_prefix}.txt', 'w') as vcf_ids_out:
            
            vcf_ids_out.write('\n'.join(self.vcf_ids))
        
        # Save variants_info
        self.variants_info.to_csv(f'{file_prefix}.tsv', sep='\t', index=False)
        
        # Save variants_matrix
        if full_matrix:
            
            np.savetxt(f'{file_prefix}.mtx', self.variants_matrix.todense(), fmt='%i', delimiter='\t')
        
        else: # Save as npz
            
            save_npz(f'{file_prefix}.npz', self.variants_matrix, compressed=True)
    
    ### ------------------------------------ ###
    
    def load_data(self, file_prefix='./vcf_data', full_matrix=False):
        
        # Read vcf_ids
        self.vcf_ids = open(f'{file_prefix}.txt', 'r').read().split('\n')
        
        # Read variants_info
        self.variants_info = pd.read_csv(f'{file_prefix}.tsv', sep='\t')
        
        # Read variants_matrix
        if full_matrix:
            
            self.variants_matrix = csc_matrix(np.loadtxt(f'{file_prefix}.mtx', dtype=int, delimiter='\t'))
        
        else: # Save as npz
            
            self.variants_matrix = load_npz(f'{file_prefix}.npz')
    
    ### ------------------------------------ ###
    ### BUILD MSA                            ###
    ### ------------------------------------ ###

    def build_phylip_msa(self, var_only=True, snp_only=False, out_name='msa.phy'):

        # Builds a Multiple Sequence Alignment using information from the vcf files
        # Note that if snp_only=True, then the +ASC option can be used for IQ-TREE
        
        # Build reference genenome using info from vcf
        chromosomes_structures, var_info_adjusted = self.build_reference_genome(self.variants_info, var_only=var_only, snp_only=snp_only)
        
        # Phylip header info
        aln_n, aln_len = self.variants_matrix.shape[1], sum([sum(bases_to_keep) for _,bases_to_keep in chromosomes_structures.values()])

        # Create MSA file
        print('Creating MSA file')
        with open(out_name, 'w') as msa_out:
            
            # Add header
            msa_out.write(f'{aln_n}\t{aln_len}\n')
            
            # Add samples
            for i in range(self.variants_matrix.shape[1]):
                
                if (i + 1) % 100 == 0:
                    
                    print(f'Processed {i + 1} / {aln_n} samples')
                
                # Get sample name
                sample_name = self.vcf_ids[i]
                
                # Extract sample's variants
                # N.B. Variants in vars_info are not sorted by position, but the ones in var_info_adjusted are!
                sample_vars_ids = self.variants_info.loc[np.ravel((self.variants_matrix[:, i] == 1).todense()), 'id'].values
                sample_vars = var_info_adjusted.loc[var_info_adjusted.id.isin(sample_vars_ids),]
                
                # Init sample's chromosomes
                sample_chroms = {key : value[0].copy() for key, value in chromosomes_structures.items()}
                
                # Updating sample's chromosomes
                for n,var in sample_vars.iterrows():
                    
                    var_id, var_type, var_chr, pos, ref, alt, _ = var
                    
                    if var_type == 'substitution':
                        
                        sample_chroms[var_chr][pos : pos + len(alt)] = list(alt)
                    
                    elif var_type == 'deletion':
                        
                        if '-' in sample_chroms[var_chr][pos : pos + len(ref)]: # Deletion overlaps with '-' areas in reference, fixing deletion size
                        
                            # Moving along the reference sequence, skipping over '-'
                            deletion_length = len(ref) - 1
                            count = 1
                            deletion_end = pos + 1
                            while count <= deletion_length:
                                
                                if sample_chroms[var_chr][deletion_end] == ref[count]:
                                    
                                    count += 1
                                    
                                deletion_end += 1
                        
                        else:
                            
                            deletion_end = pos + len(ref)
                        
                        sample_chroms[var_chr][pos + len(alt) : deletion_end] = ['-' for _ in range(deletion_end - pos - len(alt))]
                    
                    elif var_type == 'insertion':
                        
                        sample_chroms[var_chr][pos + len(ref) : pos + len(alt)] = list(alt)[len(ref):]
                    
                    else:
                        
                        pass
                
                # Remove unvariant bases if specified
                if var_only:
                    
                    for chrom in sample_chroms.keys():
                        
                        sample_chroms[chrom] = sample_chroms[chrom][chromosomes_structures[chrom][1]]
                
                # Join chromosomes into one string
                sample_chroms = ''.join([''.join(chrom) for chrom in sample_chroms.values()])
                
                # Write sample to output
                terminator = '\n' if i < aln_n - 1 else ''
                msa_out.write(f'{sample_name}\t{sample_chroms}{terminator}')
        
        print('All done!')

    ### ------------------------------------ ###

    @staticmethod
    def build_reference_genome(var_info, var_only=True, snp_only=False):
        
        # For each chromosome, create a reference based on variant sites (SNPs, insertions, and deletions)
        
        # STEPS:
        # - An chromosome composed of an array of empty strings is created
        # - Using the ref field from SNPs and deletions, info is added
        # - Where insertions occur, spacers (i.e. '-') are inserted, as many as the largest insertion at the reference base
        # - Coordinates of variants are adjusted to account for insertions
        # - If desired (var_only=True, default), only variant sites are kept and variant coordinates adjusted accordingly
        
        if snp_only:
            
            var_info = var_info.loc[var_info.type == 'substitution',]
        
        var_info_adjusted = []
        
        chromosomes_structures = {}
        for chrom in set(var_info.chr):
            
            print(f'Building reference for {chrom}')
            
            # Extract variants for chromosome
            chrom_vars = var_info.loc[var_info.chr == chrom].copy()
            
            # Sort by position
            chrom_vars.sort_values(by='pos', ascending=True, inplace=True)
            
            # Find length of chromosome (until we have variant info)
            chrom_length = max(chrom_vars.pos)
            
            # Init array of N empty strings, where N is the length of the chromosome (kinda)
            # N.B. Neither elegant nor memory efficient, but even a crappy laptop can deal with this
            chrom_struct = np.repeat('', chrom_length)
            
            # "Paint" the chromosome with info from SNPs and indels
            print('Painting SNPs and indels refs')
            for n, var in chrom_vars.loc[chrom_vars.type != 'insertion', ].iterrows():
                
                var_id, _, _, pos, ref, alt, _ = var
                
                if pos + len(ref) > chrom_length: # Elongate the chromosome due to ref info from deletions
                
                    chrom_struct = np.concatenate([chrom_struct, ['' for _ in range(pos + len(ref) - chrom_length)]])
                    chrom_length = chrom_struct.shape[0]
                
                chrom_struct[pos : pos + len(ref)] = list(ref)
            
            # Add spacers for insertions
            print('Painting insertions spacers')
            insertion_sites = chrom_vars.loc[chrom_vars.type == 'insertion', ].copy()
            insertion_pos = np.unique(insertion_sites.pos)[::-1] # Reversed so that you only affect downstream positions
            insertion_len = np.array([])
            for pos in insertion_pos: # Invert order of pos to make it easier for indexing
                
                pos_insertion = insertion_sites.loc[insertion_sites.pos == pos,] # Subset of insertion at site pos
                
                max_insertion = max([len(alt) for alt in pos_insertion.alt]) - len(ref) # Longest insertions
                insertion_len = np.append(insertion_len, max_insertion)
                
                best_ref = max([ref for ref in pos_insertion.ref], key=lambda r: len(r)) # Longest reference (should all be length 1, just being careful)
                
                chrom_struct[pos : pos + len(best_ref)] = list(best_ref) # Updating chromosome_struct with reference info from variants
                
                chrom_struct = np.concatenate([chrom_struct[:pos + len(best_ref)], ['-' for _ in range(max_insertion)], chrom_struct[pos + len(best_ref):]]) # Adding spacers for insertion
            
            # Adjust variant start positions after adding insertions
            print('Adjusting variant coordinates')
            for i_pos, i_len in zip(insertion_pos, insertion_len):
                
                chrom_vars.loc[chrom_vars.pos > i_pos, 'pos'] += int(i_len)
            
            # Trim costant sites (e.g. to make IQ-TREE faster by considering only variant sites (+ASC model))
            if var_only:
                
                print('Removing invariant sites')
                
                # Update chrom_vars
                variant_sites = np.where(chrom_struct != '')[0].tolist() # All positions that vary, be they indels, insertions, SNPs
                old_pos = np.unique(chrom_vars.pos)
                new_pos = np.where(np.isin(variant_sites, old_pos, assume_unique=True))[0]
                var_position_mapping = {o : n for o,n in zip(old_pos, new_pos)}
                chrom_vars.pos = [var_position_mapping[pos] for pos in chrom_vars.pos]
                
                # Update chrom_struct
                chrom_struct = chrom_struct[np.where(chrom_struct != '')]
                
                # Some positions that do not vary will remain due to the way insertions/deletions are marked (i.e. the first base is used for anchoring)
                # So, all positions that DO vary will here be marked for keeping
                good_pos = []
                for n, var in chrom_vars.iterrows():
                    
                    var_id, var_type, _, pos, ref, alt, _ = var
                
                    if var_type == 'substitution':
                        
                        good_pos.append(pos)
                    
                    elif var_type == 'insertion':
                        
                        good_pos.extend(list(range(pos + 1, pos + len(alt), 1)))
                    
                    elif var_type == 'deletion':
                        
                        if '-' in chrom_struct[pos : pos + len(ref)]: # Deletion overlaps with '-' areas in reference, fixing deletion size
                        
                            # Moving along the reference sequence, skipping over '-'
                            deletion_length = len(ref) - 1
                            count = 1
                            deletion_end = pos + 1
                            while count <= deletion_length:
                                
                                if chrom_struct[deletion_end] == ref[count]:
                                    
                                    count += 1
                                    
                                deletion_end += 1
                        
                        else:
                            
                            deletion_end = pos + len(ref)
                        
                        good_pos.extend(list(range(pos + 1, deletion_end, 1)))
                    
                    else:
                        
                        pass
                
                good_pos = np.unique(good_pos)
                good_paint = np.repeat(False, len(chrom_struct))
                good_paint[good_pos] = True
            
            else:
                
                good_paint = np.repeat(True, len(chrom_struct))
            
            # Convert '' to 'N' (only works for var_only=False)
            chrom_struct[chrom_struct == ''] = 'N'
            
            # Storing info
            chromosomes_structures[chrom] = [chrom_struct, good_paint]
            
            var_info_adjusted.append(chrom_vars)
        
        # Cancatenate chrom_vars
        var_info_adjusted = pd.concat(var_info_adjusted)
            
        return chromosomes_structures, var_info_adjusted

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from scipy.sparse import csc_matrix, load_npz, save_npz
