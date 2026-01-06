function initCustomQueryBuilder(form_vals) {
  const inorganic = form_vals.inorganic;
  const organic = form_vals.organic;
  const cofactors = organic.concat(inorganic);
  const evidence = form_vals.evidence;
  const tailoring = form_vals.tailoring;

  $('#builder').queryBuilder({
    plugins: ['bt-tooltip-errors'],
    allow_groups: false,
    allow_empty: false,
    operators: [
      'equal',
      'not_equal',
      'contains',
      'not_contains',
      'is_null',
      'is_not_null'
    ],
    filters: [
      {
        id: 'orcids',
        label: 'Entry ORCID',
        type: 'string',
        placeholder: 'e.g. 0000-0001-6534-6609',
        validation: {
          format: /^[0-9]/i
        },
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'references',
        label: 'Entry Reference DOI',
        type: 'string',
        placeholder: 'e.g. 10.1016/j.chembiol.2020.11.009',
        validation: {
          format: /^10\.\S+$/i
        },
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'enzyme.name',
        label: 'Enzyme Name',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.enzyme_description',
        label: 'Enzyme Description',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.uniprot_id',
        label: 'Enzyme UniProt ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.genpept_id',
        label: 'Enzyme Genpept ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.mibig_id',
        label: 'Enzyme MIBiG ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.wikidata_id',
        label: 'Enzyme Wikidata QID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.has_auxenzymes',
        label: 'Enzyme, Auxilliary',
        type: 'boolean',
        input: 'select',
        operators: ['equal', 'not_equal'],
        values: [
          true,
          false
        ]
      },
      {
        id: 'enzyme.cofactors',
        label: 'Enzyme Cofactor',
        type: 'string',
        input: 'select',
        values: cofactors,
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'enzyme.organism_id',
        label: 'Taxonomy Organism',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.domain_id',
        label: 'Taxonomy Domain',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.kingdom_id',
        label: 'Taxonomy Kingdom',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.phylum_id',
        label: 'Taxonomy Phylum',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.class_id',
        label: 'Taxonomy Class',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.order_id',
        label: 'Taxonomy Order',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.family_id',
        label: 'Taxonomy Family',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'reactions.description',
        label: 'Reaction Description',
        type: 'string',
        operators: ['contains', 'equal']
      },
      {
        id: 'reactions.reaction_smarts',
        label: 'Reaction SMARTS',
        type: 'string',
        operators: ['contains', 'equal']
      },
      {
        id: 'evidences',
        label: 'Reaction Experimental Evidence',
        type: 'string',
        input: 'select',
        values: evidence,
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'tailoring',
        label: 'Reaction Tailoring Term',
        type: 'string',
        input: 'select',
        values: tailoring,
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'reactions.rhea_id',
        label: 'Reaction Rhea ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'reactions.ec_id',
        label: 'Reaction EC Number',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'reactions.example_reactions.smiles_substrate',
        label: 'Reaction Substrate SMILES',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal',]
      },
      {
        id: 'reactions.example_reactions.products.smiles_product',
        label: 'Reaction Product SMILES',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal',]
      },
      {
        id: 'reactions.example_reactions.is_intermediate',
        label: 'Reaction, Intermediate',
        type: 'boolean',
        input: 'select',
        operators: ['equal', 'not_equal'],
        values: [
          true,
          false
        ]
      },
    ]
  });

  $('form').on('submit', function (e) {
    const rules = $('#builder').queryBuilder('getRules');
    $('#rules-input').val(JSON.stringify(rules));
  });
};
