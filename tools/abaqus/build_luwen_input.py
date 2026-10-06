#!/usr/bin/env python3
"""Assemble the original wall mechanisms; generated input still requires MOOSE precheck."""
import argparse
import collections
import hashlib
import json
from pathlib import Path
import numpy as np
import abaqus2exodus as mesh
from prepare_luwen_wall import SOURCE_SHA


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def assemble(source, stage):
    manifest = json.loads((stage / 'preparation-manifest.json').read_text())
    assert sha(source) == manifest['source_sha256'] == SOURCE_SHA
    for item in manifest['files']:
        assert sha(stage / item['path']) == item['sha256'], item['path']
    if (stage / 'wall.i').exists():
        raise ValueError('Do not overwrite an assembled input snapshot')
    model = mesh.parse_inp(source)
    gm, blocks, types, metadata, nodesets, sidesets = mesh.build_global_mesh(
        model, 1e-9, preserve_instance_nodes=True)
    gm.coords = [tuple(v * 1e-3 for v in point) for point in gm.coords]
    geometry = json.loads((stage / 'geometry-constraint-map.json').read_text())
    materials = json.loads((stage / 'materials-si.json').read_text())
    contact = json.loads((stage / 'contact-source-and-map.json').read_text())
    solid = [b for b in blocks if types[b] == 'HEX8']
    aac = [b for b in solid if metadata[b][1] == 'aac201-lunwen']
    c40 = [b for b in solid if metadata[b][1] == 'C40']
    steel = [b for b in blocks if types[b] == 'TRUSS']
    faces = []
    for card in contact['original_cards']:
        if card['keyword'] == 'contact pair':
            faces.append([item['exodus_sideset'] for item in card['mapped_surfaces']])
    # Explicit geometry-based friction stick scale. Preserve original dimensionless
    # slip tolerance; this length definition needs Abaqus runtime/reference audit.
    lengths = []
    for names in faces:
        for name in names:
            for elset, face in model.surfaces[next(k for k in model.surfaces if mesh.sanitize(k)==name)]:
                for instance, ids in model.asm_elsets[elset].items():
                    part = model.parts[mesh.model_inst_part(model, instance)]
                    for eid in ids:
                        pts = np.array([gm.coords[gm.node_map[(instance, part.elems[eid][i-1])]-1]
                                        for i in mesh.C3D8_FACES[face]])
                        lengths.append(float(np.mean(np.linalg.norm(pts-np.roll(pts,1,axis=0),axis=1))))
    characteristic_length = float(np.mean(lengths))
    slip = characteristic_length * contact['si_parameters']['slip_tolerance']
    sections = collections.OrderedDict()
    def add(section, name, **params):
        sections.setdefault(section, []).append((name, params))
    def names(values):
        return "'" + ' '.join(map(str, values)) + "'"

    add('Mesh', 'source', type='FileMeshGenerator', file='wall_mesh.e', allow_renumbering='false')
    previous = 'source'
    add('Mesh', 'rp_surface', type='LowerDBlockFromSidesetGenerator', input=previous,
        sidesets='_PickedSurf535', new_block_id=500, new_block_name='rp_surface')
    previous = 'rp_surface'
    for pair, (secondary, primary) in enumerate(faces, 1):
        for side, surface in (('secondary', secondary), ('primary', primary)):
            current = f'interface_{pair}_{side}'
            add('Mesh', current, type='LowerDBlockFromSidesetGenerator', input=previous,
                sidesets=surface, new_block_id=500+2*pair+(side=='primary'), new_block_name=current)
            previous = current
    for axis in 'xyz':
        add('Variables', 'disp_'+axis, family='LAGRANGE', order='FIRST')
    # The RP lower-dimensional surface carries Lagrange multipliers only.
    # Register an unused zero property to satisfy FEProblem material coverage;
    # no mechanics object consumes it and it contributes no stiffness.
    add('Materials', 'rp_surface_registry', type='GenericConstantMaterial', block=500,
        prop_names='rp_surface_unused_zero', prop_values=0)
    add('Materials', 'uniform_strain', type='ComputeUniformHex8IncrementalStrain', block=names(solid),
        displacements="'disp_x disp_y disp_z'", volumetric_locking_correction='false')
    for material, material_blocks in (('C40', c40), ('aac201-lunwen', aac)):
        values = materials[material]
        add('Materials', 'elasticity_'+material.replace('-','_'), type='ComputeIsotropicElasticityTensor',
            block=names(material_blocks), youngs_modulus=values['youngs_modulus_pa'], poissons_ratio=values['poissons_ratio'])
        if material == 'C40':
            add('Materials', 'c40_stress', type='ComputeFiniteStrainElasticStress', block=names(material_blocks))
        else:
            add('Materials', 'aac_stress', type='ComputeMultipleInelasticStress', block=names(material_blocks),
                inelastic_models='aac_cdp', perform_finite_strain_rotations='false')
        for axis, component in zip('xyz', range(3)):
            add('Kernels', material.replace('-', '_')+'_'+axis, type='UniformHex8StressDivergence',
                variable='disp_'+axis, component=component, block=names(material_blocks),
                displacements="'disp_x disp_y disp_z'", volumetric_locking_correction='false',
                initial_shear_modulus=values['youngs_modulus_pa']/(2*(1+values['poissons_ratio'])))
    E = materials['gangjin']['youngs_modulus_pa']
    area_by_block = {}
    for block in steel:
        origin = gm.elem_origin[blocks[block][0][0]]
        part = model.parts[mesh.model_inst_part(model, origin[0])]
        values = {section[2] for section in part.sections}
        if len(values) != 1 or None in values:
            raise ValueError('Original steel area is unresolved')
        area = values.pop()*1e-6
        area_by_block[block] = area
        add('Materials', block+'_plastic', type='SmallStrainPlasticTruss', block=block,
            displacements="'disp_x disp_y disp_z'", youngs_modulus=E,
            yield_stress=300e6, ultimate_stress=360e6, plastic_strain_knot=0.034)
        for axis, component in zip('xyz', range(3)):
            add('Kernels', block+'_'+axis, type='StressDivergenceTensorsTruss',
                block=block, variable='disp_'+axis, displacements="'disp_x disp_y disp_z'",
                component=component, area=area, use_displaced_mesh='false')
    for axis in 'xyz':
        add('BCs', 'bottom_'+axis, type='DirichletBC', variable='disp_'+axis,
            boundary='Set_3__CC_diliang_1', value=0)
    for pair, (secondary, primary) in enumerate(faces, 1):
        common = dict(secondary_boundary=secondary, primary_boundary=primary,
                      secondary_subdomain=f'interface_{pair}_secondary',
                      primary_subdomain=f'interface_{pair}_primary',
                      use_displaced_mesh='false', interpolate_normals='false',
                      correct_edge_dropping='true', ghost_point_neighbors='true')
        pressure = f'contact_pressure_{pair}'
        uo = f'contact_{pair}'
        add('Variables', pressure, family='LAGRANGE', order='FIRST', block=f'interface_{pair}_secondary')
        add('UserObjects', uo, type='SmallSlidingCohesiveContact', **common,
            displacements="'disp_x disp_y disp_z'", normal_pressure=pressure,
            stiffness="'2.559e12 1.093e12 1.093e12'", strength="'5.63e6 5.63e6 5.63e6'",
            failure_increment=0.000241, viscosity=0.0005, friction_coefficient=0.6,
            critical_elastic_slip=slip,expected_secondary_nodes=len(nodesets['SURF_'+secondary]))
        add('Constraints', f'hard_{pair}', type='SmallSlidingHardContact', **common,
            variable=pressure, secondary_variable='disp_x', primary_variable='disp_x',
            contact=uo, c=2.559e12, segment_quadrature='FOURTH')
        for axis, component in zip('xyz', range(3)):
            add('Constraints', f'traction_{pair}_{axis}', type='SmallSlidingContactTraction', **common,
                secondary_variable='disp_'+axis, primary_variable='disp_'+axis,
                contact=uo, component=component, segment_quadrature='FOURTH')
        for quantity in ('damage', 'gap', 'pressure', 'jump_n', 'jump_s', 'jump_t',
                         'traction_n', 'traction_s', 'traction_t'):
            variable = f'interface_{pair}_{quantity}'
            add('AuxVariables', variable, family='LAGRANGE', order='FIRST', block=f'interface_{pair}_secondary')
            add('AuxKernels', variable, type='SmallSlidingContactAux', variable=variable,
                boundary=secondary, contact=uo, quantity=quantity, execute_on="'timestep_end'")
    for variable, prop in (('DamageC','DamageC'),('DamageT','DamageT'),('cdp_kappa_c','cdp_kappa_c'),('cdp_kappa_t','cdp_kappa_t')):
        add('AuxVariables', variable, family='MONOMIAL', order='CONSTANT', block=names(aac))
        add('AuxKernels', variable, type='MaterialRealAux', variable=variable, property=prop,
            block=names(aac), execute_on="'initial timestep_end'")
    # Preserve tensor components in their actual definitions. In particular,
    # small total strain is not Abaqus LE, and the two CDP plastic tensors are
    # not relabeled PE/PEEQ without an independent definition check.
    tensor_outputs = (('S', 'stress', solid), ('E', 'total_strain', solid),
                      ('cdp_backbone_p', 'cdp_backbone_plastic_strain', aac),
                      ('cdp_viscous_p', 'cdp_viscous_plastic_strain', aac))
    for prefix, prop, output_blocks in tensor_outputs:
        for i, j in ((0,0),(1,1),(2,2),(0,1),(0,2),(1,2)):
            variable = f'{prefix}{i+1}{j+1}'
            add('AuxVariables', variable, family='MONOMIAL', order='CONSTANT', block=names(output_blocks))
            add('AuxKernels', variable, type='RankTwoAux', variable=variable, rank_two_tensor=prop,
                index_i=i, index_j=j, block=names(output_blocks), execute_on="'initial timestep_end'")
    for variable, prop in (('steel_axial_stress','axial_stress'),('steel_peeq','equivalent_plastic_strain'),
                           ('steel_total_strain','total_stretch'),('steel_elastic_strain','elastic_stretch'),
                           ('steel_plastic_strain','plastic_stretch')):
        add('AuxVariables', variable, family='MONOMIAL', order='CONSTANT', block=names(steel))
        add('AuxKernels', variable, type='MaterialRealAux', variable=variable, property=prop,
            block=names(steel), execute_on="'initial timestep_end'")
    add('Postprocessors','rp_traction_integral_y',type='ElementIntegralVariablePostprocessor',variable='rp_traction_y',block=500)
    add('Postprocessors','RP_RF2_N',type='ParsedPostprocessor',pp_names='rp_traction_integral_y',expression="'-rp_traction_integral_y'")
    add('Postprocessors','RP_U2_m',type='FunctionValuePostprocessor',function='rp_y_loading')
    for variable in ('DamageC','DamageT'):
        add('Postprocessors',variable+'_max',type='ElementExtremeValue',variable=variable,block=names(aac),value_type='max')
        add('Postprocessors',variable+'_aac_volume_mean',type='ElementAverageValue',variable=variable,block=names(aac))
    for material_name, output_blocks in (('aac',aac),('c40',c40)):
        for prefix in ('S','E'):
            for component in ('11','22','33','12','13','23'):
                variable = prefix+component
                add('Postprocessors',f'{variable}_{material_name}_volume_mean',type='ElementAverageValue',
                    variable=variable,block=names(output_blocks))
    add('Postprocessors','dt',type='TimestepSize')
    add('Postprocessors','newton_iterations',type='NumNonlinearIterations')
    add('Preconditioning','smp',type='SMP',full='true')
    add('Executioner/TimeStepper','',type='IterationAdaptiveDT',dt=0.01,optimal_iterations=8,
        iteration_window=2,growth_factor=1.5,cutback_factor=0.5)
    add('Outputs','csv',type='CSV',execute_on="'initial timestep_end'")
    add('Outputs','exodus',type='Exodus',execute_on="'initial timestep_end'",time_step_interval=1)
    sections['Outputs'][-1][1].update(sync_times=names([format(i/100,'.2f') for i in range(101)]),sync_only='false')
    # Main sections merge with the source-preserving AAC/embedded/RP fragments.
    lines=['# Original luwen-wxz wall; precheck and first calculation remain pending.',
           '# SI m-N-s-kg-Pa; no added loads, contacts, core columns or artificial ties.',
           '!include aac-material-fragment.i', '!include embedded-constraint-fragment.i', '']
    for section, objects in sections.items():
        lines.append('['+section+']')
        if section=='Mesh':
            lines+=['  allow_renumbering = false', '  parallel_type = REPLICATED']
        for name, params in objects:
            if name:lines.append('  ['+name+']')
            indent='    ' if name else '  '
            for key,value in params.items():lines.append(indent+key+' = '+str(value))
            if name:lines.append('  []')
        lines+=['[]','']
        if section == 'Variables':
            # Field-variable registration must precede the RP scalar for the
            # locked framework's scaling-factor vector layout.
            lines += ['!include rp-constraint-fragment.i', '']
    custom=[geometry['blocks'][b]['id'] for b in solid]
    # Embedded and nonlocal mortar couplings extend the mesh adjacency graph.
    # Upstream embedded-constraint tests use hash assembly for this pattern.
    lines+=['[Problem]', '  type = FEProblem',
            '  use_hash_table_matrix_assembly = true', '[]', '', '[Executioner/Quadrature]',
            '  type = GAUSS', '  order = SECOND',
            '  custom_blocks = '+names(custom),
            '  custom_orders = '+names(['FIRST']*len(custom)), '[]','',
            '[Executioner]','  type = Transient','  solve_type = NEWTON','  line_search = basic',
            '  start_time = 0','  end_time = 1','  dtmin = 1e-8','  dtmax = 1',
            '  num_steps = 10000','  nl_max_its = 60','  nl_rel_tol = 1e-8','  nl_abs_tol = 1e-7',
            '  automatic_scaling = true','  compute_scaling_once = false',
            "  petsc_options_iname = '-pc_type -pc_factor_mat_solver_type'",
            "  petsc_options_value = 'lu mumps'",'[]','',
            '[Outputs]','  file_base = wall_out','  print_linear_residuals = false', '[]','']
    (stage/'wall.i').write_text('\n'.join(lines))
    record={'status':'assembled-not-prechecked','source_sha256':SOURCE_SHA,
            'assembly_generator_sha256':sha(__file__),'steel_area_m2_by_block':area_by_block,
            'friction_characteristic_length_m':characteristic_length,
            'friction_length_definition':'mean of quadrilateral contact-face mean edge lengths; Abaqus internal characteristic-length equivalence unverified',
            'critical_elastic_slip_m':slip,'contact_pairs':len(faces),'solid_material_points':14900,
            'original_physics_time_s':1,'rp_target_U2_m':-0.020,
            'limitations':['Abaqus reference absent; no equivalence claim',
                           'custom mortar discretization requires force-transfer/coverage runtime checks',
                           'C3D8R hourglass normalization requires independent reference audit'],
            'files':[{'name':p.name,'sha256':sha(p),'size':p.stat().st_size}
                     for p in sorted(stage.iterdir()) if p.is_file()]}
    (stage/'full-input-manifest.json').write_text(json.dumps(record,ensure_ascii=False,indent=2)+'\n')
    return record

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inp',required=True,type=Path)
    parser.add_argument('--stage',required=True,type=Path)
    args=parser.parse_args()
    result=assemble(args.inp.resolve(),args.stage.resolve())
    print(json.dumps({k:result[k] for k in ('status','contact_pairs','critical_elastic_slip_m')}))
