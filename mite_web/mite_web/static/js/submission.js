function removeEntry(button, class_ref) {
    const targetDiv = button.closest(class_ref);
    if (targetDiv) {
        targetDiv.remove();
    }
}

function addEnzymeRef() {
    const container = document.getElementById('enzyme-ref');
    const entryHtml = `
        <div class="enzyme_ref">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="enzyme_ref[]" id="enzyme_ref[]" class="form-control" aria-describedby="Enzyme Reference DOI Help"  value="" required>
                        <label for="enzyme_ref[]" class="form-label">Enzyme Reference DOI</label>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.enzyme_ref')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}


function addReactionRef(index) {
    const container = document.getElementById(index);
    const entryHtml = `
        <div class="reaction_ref">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${index}" id="${index}" class="form-control" aria-describedby="Reaction Reference DOI Help"  value="" required>
                        <label for="${index}" class="form-label">Reaction Reference DOI</label>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.reaction_ref')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function addReactionProd(index) {
    const container = document.getElementById(index);
    const entryHtml = `
        <div class="reaction_products">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${index}" id="${index}" class="form-control" aria-describedby="${index}Help"  value="" required>
                        <label for="${index}" class="form-label">Product SMILES</label>
                        <div id="${index}Help" class="form-text">A single reaction product SMILES string. Dot notation is not permitted</div>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.reaction_products')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function addForbiddenProd(index) {
    const container = document.getElementById(index);
    const entryHtml = `
        <div class="forbidden_products">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${index}" id="${index}" class="form-control" aria-describedby="${index}Help"  value="" required>
                        <label for="${index}" class="form-label">Forbidden Product SMILES</label>
                        <div id="${index}Help" class="form-text">A single forbidden reaction product SMILES string. Dot notation is not permitted</div>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.forbidden_products')">Remove</button>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function addAuxEnzyme() {
    const container = document.getElementById('aux-enzyme');
    const index = container.children.length;
    const entryHtml = `
        <div class="aux_enzyme">
            <div class="card card-body">
                <div class="row mb-2">
                    <div class="col">
                        <h6 class="fw-semibold lh-2">Auxiliary Enzyme</h6>
                    </div>
                    <div class="col-auto mx-auto">
                        <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.aux_enzyme')">Remove</button>
                    </div>
                </div>
                <div class="row">
                    <div class="col">
                         <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]name" id="auxenzyme[${index}]name" class="form-control" aria-describedby="AuxEnzymeNameHelp"  value="" required>
                            <label for="auxenzyme[${index}]name" class="form-label">Auxiliary Enyzme Name</label>
                            <div id="AuxEnzymeNameHelp" class="form-text">The common enzyme name (e.g. McjC)</div>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]description" id="auxenzyme[${index}]description" class="form-control" aria-describedby="AuxEnzymeDescriptionHelp"  value="" required>
                            <label for="auxenzyme[${index}]description" class="form-label">Auxiliary Enzyme Description</label>
                            <div id="AuxEnzymeDescriptionHelp" class="form-text">A brief description of the enzyme function</div>
                        </div>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]uniprot" id="auxenzyme[${index}]uniprot" class="form-control" aria-describedby="AuxEnzymeUniprotHelp"  value="">
                            <label for="auxenzyme[${index}]uniprot" class="form-label">UniProt/UniParc ID</label>
                            <div id="AuxEnzymeUniprotHelp" class="form-text">The UniProtKB or UniParc ID</div>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]genpept" id="auxenzyme[${index}]genpept" class="form-control" aria-describedby="AuxEnzymeGenpeptHelp"  value="">
                            <label for="auxenzyme[${index}]genpept" class="form-label">GenPept ID</label>
                            <div id="AuxEnzymeGenpeptHelp" class="form-text">The NCBI GenPept ID</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}


function addKnownReaction(knownreaction) {
    const container = document.getElementById(knownreaction);
    const index = container.children.length;
    const entryHtml = `
        <div class="knownreaction">
            <div class="card card-body">
                <div class="row mb-2">
                    <div class="col">
                        <h6 class="lh-2">Known Reaction</h6>
                    </div>
                    <div class="col-auto mx-auto">
                        <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.knownreaction')">Remove</button>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="${knownreaction}substrate" id="${knownreaction}substrate" class="form-control" aria-describedby="KnownreactionSubstrateHelp"  value="" required>
                            <label for="${knownreaction}substrate" class="form-label">Substrate SMILES</label>
                            <div id="KnownreactionSubstrateHelp" class="form-text">The reaction substrate SMILES string. Multiple substrates must be specified with dot notation ('substrate1.substrate2')</div>
                        </div>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="${knownreaction}description" id="${knownreaction}description" class="form-control" aria-describedby="KnownreactionDescriptionHelp"  value="">
                            <label for="${knownreaction}description" class="form-label">Description</label>
                            <div id="KnownreactionDescriptionHelp" class="form-text">Description of the reaction example</div>
                        </div>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col-auto">
                        <div class="form-check-inline">
                            <input class="form-check-input" value="True" type="radio" name='${knownreaction}intermediate' id='${knownreaction}intermediate1'>
                            <label class="form-check-label" for='${knownreaction}intermediate1'>True</label>
                        </div>
                        <div class="form-check-inline">
                            <input class="form-check-input" value="False" type="radio" name='${knownreaction}intermediate' id='${knownreaction}intermediate2' checked>
                            <label class="form-check-label" for='${knownreaction}intermediate2'>False</label>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-text">Intermediate product?</div>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div id='${knownreaction}products[]'>
                        </div>
                        <button type="button" class="btn btn-secondary my-2" onclick="addReactionProd('${knownreaction}products[]')">Add Product SMILES</button>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div id='{{ "reaction[" ~ reaction_outer_idx ~ "]knownreaction[" ~ reaction_example_idx ~ "]forbiddenproducts[]" }}'>
                        {% if knownreaction.get("forbidden_products") %}
                        {% for product in knownreaction.forbidden_products %}
                        <div class="forbidden_products">
                            <div class="row d-flex align-items-center mb-1">
                                <div class="col">
                                    {{ macro.render_text(title="Forbidden Product SMILES", id="reaction[" ~ reaction_outer_idx ~ "]knownreaction[" ~ reaction_example_idx ~ "]forbiddenproducts[]", type="text", value=product, required=true, help="A single forbidden reaction product SMILES string. Dot notation is not permitted") }}
                                </div>
                                <div class="col-auto mx-auto">
                                    <button type="button" class="btn btn-danger" onclick="removeEntry(this, '.forbidden_products')">Remove</button>
                                </div>
                            </div>
                        </div>
                        {% endfor %}
                        {% endif %}
                        </div>
                        <button type="button" class="btn btn-secondary my-2" onclick="addForbiddenProd('{{ 'reaction[' ~ reaction_outer_idx ~ ']knownreaction[' ~ reaction_example_idx ~ ']forbiddenproducts[]' }}')">Add Forbidden Product SMILES</button>
                    </div>
                </div>
            </div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}
