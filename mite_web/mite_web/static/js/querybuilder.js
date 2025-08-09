function initCustomQueryBuilder(form_vals) {
  const inorganic = form_vals.inorganic;
  const organic = form_vals.organic;
  const cofactors = organic.concat(inorganic);

  $('#builder').queryBuilder({
    plugins: ['bt-tooltip-errors'],
    allow_groups: false,
    allow_empty: true,
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
        label: 'ORCID',
        type: 'string',
        placeholder: 'e.g. 0000-0001-6534-6609',
        validation: {
          format: /^[0-9]/i
        },
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },
      {
        id: 'references',
        label: 'Literature DOI',
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
        label: 'UniProt ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.genpept_id',
        label: 'Genpept ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.mibig_id',
        label: 'MIBiG ID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.wikidata_id',
        label: 'Wikidata QID',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.has_auxenzymes',
        label: 'Auxilliary Enzymes',
        type: 'boolean',
        input: 'select',
        operators: ['equal', 'not_equal'],
        values: [
          true,
          false
        ]
      },
      {
        id: 'enzyme.organism_id',
        label: 'Organism',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.domain_id',
        label: 'Domain',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.kingdom_id',
        label: 'Kingdom',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.phylum_id',
        label: 'Phylum',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.class_id',
        label: 'Class',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.order_id',
        label: 'Order',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.family_id',
        label: 'Family',
        type: 'string',
        operators: ['contains', 'not_contains', 'equal', 'not_equal', 'is_null', 'is_not_null']
      },
      {
        id: 'enzyme.cofactors',
        label: 'Cofactor',
        type: 'string',
        input: 'select',
        values: cofactors,
        validation: {
          format: /^10\.\S+$/i
        },
        operators: ['contains', 'not_contains', 'equal', 'not_equal']
      },






    ]
  });

  $('#builder').queryBuilder('setRules', {
      condition: 'AND',
      rules: []
  });

  $('form').on('submit', function (e) {
    const rules = $('#builder').queryBuilder('getRules');
    if (!rules || !rules.valid) {
      alert("Invalid query!");
      e.preventDefault(); // stop form from submitting
      return;
    }
    $('#rules-input').val(JSON.stringify(rules));
  });
};
