// Validates a simple addition on form submission, block submission if invalid
function validateSumInput(event, x, y) {
    const userSum = parseInt(document.getElementById('usersum').value, 10);
    const correctSum = +x + +y;

    const sumInputForm = document.getElementById('usersum');

    if (userSum !== correctSum) {
        sumInputForm.style.borderColor = 'red';
        document.getElementById('usersum-error-message').textContent = "Incorrect sum! Please check your answer.";
        event.preventDefault();
    } else {
        sumInputForm.style.borderColor = '';
        document.getElementById('usersum-error-message').textContent = "";
    }
}



// Removes a Div of a specified class closest to the button; used by all elements that can be added/removed by user-input
function removeField(button, class_ref) {
    const targetDiv = button.closest(class_ref);
    if (targetDiv) {
        targetDiv.remove();
    }
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeName(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_name" name="enzyme_name" class="form-control" aria-describedby="EnzymeNameHelp" type="text" value='${data?.enzyme?.name ?? ""}' required>
            <label for="enzyme_name" class="form-label">Enzyme Name</label>
            <div id="EnzymeNameHelp" class="form-text">The common enzyme name (e.g. McjB)</div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeDescription(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_description" name="enzyme_description" class="form-control" aria-describedby="EnzymeDescriptionHelp" type="text" value='${data?.enzyme?.description ?? ""}' required>
            <label for="enzyme_description" class="form-label">Enzyme Description</label>
            <div id="EnzymeDescriptionHelp" class="form-text">A brief description of the enzyme function</div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeUniprot(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_uniprot" name="enzyme_uniprot" class="form-control" aria-describedby="EnzymeUniprotHelp" type="text" value='${data?.enzyme?.databaseIds?.uniprot ?? ""}'>
            <label for="enzyme_description" class="form-label">Enzyme UniProt ID</label>
            <div id="EnzymeUniprotHelp" class="form-text">The UniProtKB or UniParc ID</div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeGenpept(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_genpept" name="enzyme_genpept" class="form-control" aria-describedby="EnzymeGenpeptHelp" type="text" value='${data?.enzyme?.databaseIds?.genpept ?? ""}'>
            <label for="enzyme_genpept" class="form-label">Enzyme GenPept ID</label>
            <div id="EnzymeGenpeptHelp" class="form-text">The NCBI GenPept ID</div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate a unique field, triggered by DOMContentLoaded event
function createHtmlEnzymeMibig(container, data = {}) {
    const entryHtml = `
        <div class="form-floating">
            <input id="enzyme_mibig" name="enzyme_mibig" pattern="^BGC.+" minlength="10" maxlength="10" class="form-control" aria-describedby="EnzymeMibigHelp" type="text" value='${data?.enzyme?.databaseIds?.mibig ?? ""}'>
            <label for="enzyme_mibig" class="form-label">MIBiG ID</label>
            <div id="EnzymeMibigHelp" class="form-text">The MIBiG ID of the BGC containing the enzyme</div>
        </div>
    `;
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Populate the enzyme reference forms, triggered by DOMContentLoaded event
function populateEnzymeRefForm(container, data) {
    if (data.enzyme?.references) {
        data.enzyme.references.forEach(entry => {
            const entryHtml = createHtmlRef(entry, "enzyme_ref", "enzyme_ref[]");
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

// Populate the auxiliary Enzyme forms, triggered by DOMContentLoaded event
function populateAuxEnzymeForm(container, data) {
    if (data.enzyme?.auxiliaryEnzymes) {
        data.enzyme.auxiliaryEnzymes.forEach(entry => {
            const index = container.children.length;
            const entryHtml = createHtmlAuxEnyzme(entry, index)
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

// Populate the Reaction forms, triggered by DOMContentLoaded event
function populateReactionForm(container, data) {
    if (data.reactions) {
        data.reactions.forEach(entry => {
            const index = container.children.length;
            createHtmlReaction(entry, index)
        });
    }
}

// Triggered by button to add an empty reference form
function insertEnzymeRefForm(container) {
    const entryHtml = createHtmlRef("", "enzyme_ref", "enzyme_ref[]")
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggered by button to add an empty AuxEnzyme
function insertAuxEnzymeForm(container) {
    const index = container.children.length;
    const entryHtml = createHtmlAuxEnyzme({}, index);
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggered by button to add an empty ReactionForm
function insertReactionForm(container) {
    const index = container.children.length;
    createHtmlReaction({}, index)
}

// Triggred by button to add an empty ReactionReferenceForm
function insertReactionRefForm(index) {
    const container = document.getElementById(`reaction[${index}]ref-field`);
    const entryHtml = createHtmlRef("", "reaction_ref", `reaction[${index}]ref[]`);
    container.insertAdjacentHTML('beforeend', entryHtml);
}

// Triggred by button to add an empty KnownReactionForm
function insertKnownReactionForm(index) {
    const knownReactionElements = document.querySelectorAll('.knownreaction');
    const containerReaction = document.getElementById(`reaction[${index}]knownreaction-field`);
    const index_inner = containerReaction.children.length;
    addHtmlKnownReaction({}, index, index_inner);
}

// Creates an reference field, used by enzyme information and Reaction fields
function createHtmlRef(data, className, jsonID) {
    return `
        <div class="${className}">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${jsonID}" pattern="^doi:10.+" id="${jsonID}" class="form-control" aria-describedby="${jsonID}Help"  value="${data}" required>
                        <label for="${jsonID}" class="form-label">Reference DOI</label>
                        <div id="${jsonID}Help" class="form-text">A Digital Object Identifier (DOI), must begin with 'doi:10.'</div>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.${className}')">Remove</button>
                </div>
            </div>
        </div>
    `;
}

// Creates the auxilary enzyme field; index is used to keep track of the instance that is opened
function createHtmlAuxEnyzme(data = {}, index) {
    return `
        <div class="aux_enzyme">
            <div class="card card-body">
                <div class="row mb-2">
                    <div class="col">
                        <h6 class="fw-semibold lh-2">Auxiliary Enzyme</h6>
                    </div>
                    <div class="col-auto mx-auto">
                        <button type="button" class="btn btn-danger" onclick="removeField(this, '.aux_enzyme')">Remove</button>
                    </div>
                </div>
                <div class="row">
                    <div class="col">
                         <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]name" id="auxenzyme[${index}]name" class="form-control" aria-describedby="AuxEnzymeNameHelp"  value='${data?.name ?? ""}' required>
                            <label for="auxenzyme[${index}]name" class="form-label">Auxiliary Enyzme Name</label>
                            <div id="AuxEnzymeNameHelp" class="form-text">The common enzyme name (e.g. McjC)</div>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]description" id="auxenzyme[${index}]description" class="form-control" aria-describedby="AuxEnzymeDescriptionHelp"  value='${data?.description ?? ""}' required>
                            <label for="auxenzyme[${index}]description" class="form-label">Auxiliary Enzyme Description</label>
                            <div id="AuxEnzymeDescriptionHelp" class="form-text">A brief description of the enzyme function</div>
                        </div>
                    </div>
                </div>
                <div class="row mb-2">
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]uniprot" id="auxenzyme[${index}]uniprot" class="form-control" aria-describedby="AuxEnzymeUniprotHelp"  value='${data?.databaseIds?.uniprot ?? ""}'>
                            <label for="auxenzyme[${index}]uniprot" class="form-label">UniProt/UniParc ID</label>
                            <div id="AuxEnzymeUniprotHelp" class="form-text">The UniProtKB or UniParc ID</div>
                        </div>
                    </div>
                    <div class="col">
                        <div class="form-floating">
                            <input type="text" name="auxenzyme[${index}]genpept" id="auxenzyme[${index}]genpept" class="form-control" aria-describedby="AuxEnzymeGenpeptHelp"  value='${data?.databaseIds?.genpept ?? ""}'>
                            <label for="auxenzyme[${index}]genpept" class="form-label">GenPept ID</label>
                            <div id="AuxEnzymeGenpeptHelp" class="form-text">The NCBI GenPept ID</div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
}

// Creates the reaction field; index is used to keep track of the instance that is opened
function createHtmlReaction(data = {}, index) {
    const reactionEntry = document.createElement('div')
    reactionEntry.classList.add('reaction')
    reactionEntry.innerHTML = `
        <div class="card card-body">
            <div class="row mb-2">
                <div class="col">
                    <h5 class="fw-semibold lh-2">Reaction</h5>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.reaction')">Remove</button>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="card card-body">
                        <div class="row mb-2">
                            <div class="col">
                                <h6 class="fw-semibold lh-2">General Information</h6>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]smarts" name="reaction[${index}]smarts" class="form-control" aria-describedby="ReactionSmartsHelp" type="text" value='${data?.reactionSMARTS ?? ""}' required>
                                    <label for="reaction[${index}]smart" class="form-label">Reaction SMARTS</label>
                                    <div id="ReactionSmartsHelp" class="form-text">The reaction SMARTS depicting the reaction</div>
                                </div>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]description" name="reaction[${index}]description" class="form-control" aria-describedby="ReactionDescriptionHelp" type="text" value='${data?.description ?? ""}' required>
                                    <label for="reaction[${index}]description" class="form-label">Description</label>
                                    <div id="ReactionDescriptionHelp" class="form-text">Brief description of the reaction</div>
                                </div>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]rhea" name="reaction[${index}]rhea" class="form-control" aria-describedby="ReactionRheaHelp" type="text" value='${data?.databaseIds?.rhea ?? ""}'>
                                    <label for="reaction[${index}]rhea" class="form-label">RHEA ID</label>
                                    <div id="ReactionRheaHelp" class="form-text">The RHEA Knowledgebase ID</div>
                                </div>
                            </div>
                            <div class="col">
                                <div class="form-floating">
                                    <input id="reaction[${index}]ec" name="reaction[${index}]ec" pattern="^[0-9]+\..+" class="form-control" aria-describedby="ReactionEcHelp" type="text" value='${data?.databaseIds?.ec ?? ""}'>
                                    <label for="reaction[${index}]ec" class="form-label">EC Number</label>
                                    <div id="ReactionEcHelp" class="form-text">The EC (Enzyme Commission) Number (e.g. '1.14.19.59')</div>
                                </div>
                            </div>
                        </div>
                        <div class="row">
                            <div class="col">
                                <h6 class="lh-2">Tailoring Reaction Controlled Vocabulary</h6>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div class="form-group">
                                    <div class="card card-body">
                                        <div class="row" id="reaction[${index}]tailoring-field">
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="row">
                            <div class="col">
                                <h6 class="lh-2">Experimental Evidence Qualifiers</h6>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div class="form-group">
                                    <div class="card card-body">
                                        <div class="row" id="reaction[${index}]evidencecode-field">
                                        </div>
                                    </div>
                                </div>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <div id="reaction[${index}]ref-field"></div>
                            </div>
                        </div>
                        <div class="row mb-2">
                            <div class="col">
                                <button type="button" class="btn btn-secondary my-2" onclick="insertReactionRefForm('${index}')">Add Reaction Reference</button>
                            </div>
                        </div>
                    </div>
                </div>
                <div class="col">
                    <div class="card card-body">
                        <div class="row mb-2">
                            <div class="col">
                                <h6 class="fw-semibold lh-2">Known Reactions</h6>
                            </div>
                        </div>
                        <div id="reaction[${index}]knownreaction-field"></div>
                        <div class="row">
                            <div class="col">
                                <button type="button" class="btn btn-secondary my-2" onclick="insertKnownReactionForm(${index})">Add Known Reaction</button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;
    const container = document.getElementById('reactions');
    container.appendChild(reactionEntry);

    let tailoringArray = [
        "Acetylation",
        "Acylation",
        "Amination",
        "Biaryl bond formation",
        "Carboxylation",
        "Cyclization",
        "Deamination",
        "Decarboxylation",
        "Dehydration",
        "Dehydrogenation",
        "Demethylation",
        "Dioxygenation",
        "Epimerization",
        "Glycosylation",
        "Halogenation",
        "Heterocyclization",
        "Hydrolysis",
        "Hydroxylation",
        "Macrolactam formation",
        "Methylation",
        "Monooxygenation",
        "Oxidation",
        "Phosphorylation",
        "Prenylation",
        "Reduction",
        "Sulfation",
        "Other"
    ];
    const containerTailoring = document.getElementById(`reaction[${index}]tailoring-field`);
    for (let tailoring of tailoringArray) {
        if (data.tailoring && data.tailoring.includes(tailoring) ) {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col');
            entryHtml.innerHTML = addTickedCheckbox(`reaction[${index}]tailoring[]`, tailoring);
            containerTailoring.appendChild(entryHtml)
        } else {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col');
            entryHtml.innerHTML =  addUntickedCheckbox(`reaction[${index}]tailoring[]`, tailoring);
            containerTailoring.appendChild(entryHtml)
        }
    }

    let evidenceCodeArray = [
        "Heterologous expression",
        "Inference from genomic data and chemical structure",
        "In vitro assay",
        "Isothermal titration calorimetry",
        "Knock-out studies",
        "Site-directed mutagenesis",
        "Surface plasmon resonance"
    ]
    const containerEvidenceCode = document.getElementById(`reaction[${index}]evidencecode-field`);
    for (let evidence of evidenceCodeArray) {
        if (data.evidence?.evidenceCode && data.evidence.evidenceCode.includes(evidence) ) {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col');
            entryHtml.innerHTML = addTickedCheckbox(`reaction[${index}]evidencecode[]`, evidence);
            containerEvidenceCode.appendChild(entryHtml);
        } else {
            const entryHtml = document.createElement('div');
            entryHtml.classList.add('col');
            entryHtml.innerHTML =  addUntickedCheckbox(`reaction[${index}]evidencecode[]`, evidence);
            containerEvidenceCode.appendChild(entryHtml);
        }
    };

    const containerEvidenceRef = document.getElementById(`reaction[${index}]ref-field`);
    if (data.evidence?.references) {
        data.evidence.references.forEach(ref => {
            const entryHtml = document.createElement('div');
            entryHtml.innerHTML =  createHtmlRef(ref, "reaction_ref", `reaction[${index}]ref[]`);
            containerEvidenceRef.appendChild(entryHtml);
        });
    };

    if (data.reactions) {
        data.reactions.forEach(data_reaction => {
            const containerKR = document.getElementById(`reaction[${index}]knownreaction-field`);
            const indexKnownReaction = containerKR.children.length;
            addHtmlKnownReaction(data_reaction, index, indexKnownReaction);
        });
    };
}

// Adds a ticked checkbox with the appropriate value
function addTickedCheckbox(container_id, qualifier) {
    return `
        <div class="form-check">
            <input class="form-check-input" type="checkbox" id="${qualifier}" name='${container_id}' value="${qualifier}" checked>
            <label class="form-check-label" for="${qualifier}">${qualifier}</label>
        </div>
    `;
}

// Adds an unticked checkbox with the appropriate value
function addUntickedCheckbox(container_id, qualifier) {
    return `
        <div class="form-check">
            <input class="form-check-input" type="checkbox" id="${qualifier}" name='${container_id}' value="${qualifier}">
            <label class="form-check-label" for="${qualifier}">${qualifier}</label>
        </div>
    `;
}

// Inserts a KnownReaction Form into an existing Reaction form; outer index tracks the instance of the Reaction, the inner index the instance of the Known reaction inside the Reaction
function addHtmlKnownReaction(data, index_outer, index_inner) {
    const knownReactionEntry = document.createElement('div');
    knownReactionEntry.classList.add('knownreaction');
    knownReactionEntry.innerHTML = `
        <div class="card card-body">
            <div class="row mb-2">
                <div class="col">
                    <h6 class="lh-2">Known Reaction</h6>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.knownreaction')">Remove</button>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="reaction[${index_outer}]knownreaction[${index_inner}]substrate" id="reaction[${index_outer}]knownreaction[${index_inner}]substrate" class="form-control" aria-describedby="KnownReactionSubstrateHelp"  value='${data?.substrate ?? ""}' required>
                        <label for="reaction[${index_outer}]knownreaction[${index_inner}]substrate" class="form-label">Substrate SMILES</label>
                        <div id="KnownReactionSubstrateHelp" class="form-text">The reaction substrate SMILES string. Multiple substrates must be specified with dot notation ('substrate1.substrate2')</div>
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="reaction[${index_outer}]knownreaction[${index_inner}]description" id="reaction[${index_outer}]knownreaction[${index_inner}]description" class="form-control" aria-describedby="KnownReactionDescriptionHelp"  value='${data?.description ?? ""}'>
                        <label for="reaction[${index_outer}]knownreaction[${index_inner}]description" class="form-label">Description</label>
                        <div id="KnownReactionSubstrateHelp" class="form-text">Description of the reaction example</div>
                    </div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col-auto">
                    <div class="form-check-inline">
                        <input class="form-check-input" value="True" type="radio" name="reaction[${index_outer}]knownreaction[${index_inner}]intermediate" id="reaction[${index_outer}]knownreaction[${index_inner}]intermediate1">
                        <label class="form-check-label" for="reaction[${index_outer}]knownreaction[${index_inner}]intermediate1">True</label>
                    </div>
                    <div class="form-check-inline">
                        <input class="form-check-input" value="False" type="radio" name="reaction[${index_outer}]knownreaction[${index_inner}]intermediate" id="reaction[${index_outer}]knownreaction[${index_inner}]intermediate2" checked>
                        <label class="form-check-label" for="reaction[${index_outer}]knownreaction[${index_inner}]intermediate2">False</label>
                    </div>
                </div>
                <div class="col">
                    <div class="form-text">Intermediate product?</div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div id="reaction[${index_outer}]knownreaction[${index_inner}]products[]"></div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <button type="button" class="btn btn-secondary my-2" onclick="insertProductForm('reaction[${index_outer}]knownreaction[${index_inner}]products[]')">Add Product</button>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <div id="reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]"></div>
                </div>
            </div>
            <div class="row mb-2">
                <div class="col">
                    <button type="button" class="btn btn-secondary my-2" onclick="insertForbiddenProductForm('reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]')">Add Forbidden Product</button>
                </div>
            </div>
        </div>
    `;
    const containerReaction = document.getElementById(`reaction[${index_outer}]knownreaction-field`);
    containerReaction.appendChild(knownReactionEntry);

    //    populates product fields based on input data
    if (data.products) {
        data.products.forEach(product => {
            insertProductForm(`reaction[${index_outer}]knownreaction[${index_inner}]products[]`, product);
        });
    };

    //    populates forbidden product fields based on input data
    if (data.forbidden_products) {
        data.forbidden_products.forEach(product => {
            insertForbiddenProductForm(`reaction[${index_outer}]knownreaction[${index_inner}]forbiddenproducts[]`, product);
        });
    };

}

// Inserts a knownProduct form inside the Known Reactions field
function insertProductForm(container_id, data = "") {
    const reactionProductEntry = document.createElement('div');
    reactionProductEntry.classList.add('reactionproducts');
    reactionProductEntry.innerHTML = `
        <div class="row d-flex align-items-center mb-1">
            <div class="col">
                <div class="form-floating">
                    <input type="text" name="${container_id}" id="${container_id}" class="form-control" aria-describedby="ReactionProductsHelp" value='${data}' required>
                    <label for="${container_id}" class="form-label">Product SMILES</label>
                    <div id="ReactionProductsHelp" class="form-text">A single reaction product SMILES string. Dot notation is not permitted</div>
                </div>
            </div>
            <div class="col-auto mx-auto">
                <button type="button" class="btn btn-danger" onclick="removeField(this, '.reactionproducts')">Remove</button>
            </div>
        </div>
    `;
    const container = document.getElementById(container_id);
    container.appendChild(reactionProductEntry);
}

// Inserts a ForbiddenProduct form inside the Known Reactions field
function insertForbiddenProductForm(container_id, data = "") {
    const forbiddenProductEntry = document.createElement('div');
    forbiddenProductEntry.classList.add('forbiddenproducts');
    forbiddenProductEntry.innerHTML = `
        <div class="row d-flex align-items-center mb-1">
            <div class="col">
                <div class="form-floating">
                    <input type="text" name="${container_id}" id="${container_id}" class="form-control" aria-describedby="ForbiddenProductsHelp" value='${data}' required>
                    <label for="${container_id}" class="form-label">Forbidden Product SMILES</label>
                    <div id="ForbiddenProductsHelp" class="form-text">A single forbidden reaction product SMILES string. Dot notation is not permitted</div>
                </div>
            </div>
            <div class="col-auto mx-auto">
                <button type="button" class="btn btn-danger" onclick="removeField(this, '.forbiddenproducts')">Remove</button>
            </div>
        </div>
    `;
    const container = document.getElementById(container_id);
    container.appendChild(forbiddenProductEntry);
}