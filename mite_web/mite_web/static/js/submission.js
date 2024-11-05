// Removal function

function removeField(button, class_ref) {
    const targetDiv = button.closest(class_ref);
    if (targetDiv) {
        targetDiv.remove();
    }
}

// DOM Scripts (populate fields using input data)

function populateTextForm(container, htmlFunction, data) {
    const entryHtml = htmlFunction(data);
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function populateEnzymeRefForm(container, data) {
    if (data.enzyme?.references) {
        data.enzyme.references.forEach(entry => {
            const entryHtml = createHtmlRef(entry, "enzyme_ref", "enzyme_ref[]");
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

function populateAuxEnzymeForm(container, data) {
    if (data.enzyme?.auxiliaryEnzymes) {
        data.enzyme.auxiliaryEnzymes.forEach(entry => {
            const index = container.children.length;
            const entryHtml = createHtmlAuxEnyzme(entry, index)
            container.insertAdjacentHTML('beforeend', entryHtml);
        });
    }
}

function populateReactionForm(container, data) {
    if (data.reactions) {
        data.reactions.forEach(entry => {
            const index = container.children.length;
            createHtmlReaction(entry, index)
        });
    }
}

// On-demand scripts (insert forms)

function insertEnzymeRefForm(container) {
    const entryHtml = createHtmlRef("", "enzyme_ref", "enzyme_ref[]")
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function insertAuxEnzymeForm(container) {
    const index = container.children.length;
    const entryHtml = createHtmlAuxEnyzme({}, index)
    container.insertAdjacentHTML('beforeend', entryHtml);
}

function insertReactionForm(container) {
    const index = container.children.length;
    createHtmlReaction({}, index)
}


// All HTML generation scripts


function createHtmlEnzymeName(data = {}) {
    return `
        <div class="form-floating">
            <input id="enzyme_name" name="enzyme_name" class="form-control" aria-describedby="EnzymeNameHelp" type="text" value='${data?.enzyme?.name ?? ""}' required>
            <label for="enzyme_name" class="form-label">Enzyme Name</label>
            <div id="EnzymeNameHelp" class="form-text">The common enzyme name (e.g. McjB)</div>
        </div>
    `;
}

function createHtmlEnzymeDescription(data = {}) {
    return `
        <div class="form-floating">
            <input id="enzyme_description" name="enzyme_description" class="form-control" aria-describedby="EnzymeDescriptionHelp" type="text" value='${data?.enzyme?.description ?? ""}' required>
            <label for="enzyme_description" class="form-label">Enzyme Description</label>
            <div id="EnzymeDescriptionHelp" class="form-text">A brief description of the enzyme function</div>
        </div>
    `;
}

function createHtmlEnzymeUniprot(data = {}) {
    return `
        <div class="form-floating">
            <input id="enzyme_uniprot" name="enzyme_uniprot" class="form-control" aria-describedby="EnzymeUniprotHelp" type="text" value='${data?.enzyme?.databaseIds?.uniprot ?? ""}'>
            <label for="enzyme_description" class="form-label">Enzyme UniProt ID</label>
            <div id="EnzymeUniprotHelp" class="form-text">The UniProtKB or UniParc ID</div>
        </div>
    `;
}

function createHtmlEnzymeGenpept(data = {}) {
    return `
        <div class="form-floating">
            <input id="enzyme_genpept" name="enzyme_genpept" class="form-control" aria-describedby="EnzymeGenpeptHelp" type="text" value='${data?.enzyme?.databaseIds?.genpept ?? ""}'>
            <label for="enzyme_genpept" class="form-label">Enzyme GenPept ID</label>
            <div id="EnzymeGenpeptHelp" class="form-text">The NCBI GenPept ID</div>
        </div>
    `;
}

function createHtmlEnzymeMibig(data = {}) {
    return `
        <div class="form-floating">
            <input id="enzyme_mibig" name="enzyme_mibig" class="form-control" aria-describedby="EnzymeMibigHelp" type="text" value='${data?.enzyme?.databaseIds?.mibig ?? ""}'>
            <label for="enzyme_mibig" class="form-label">MIBiG ID</label>
            <div id="EnzymeMibigHelp" class="form-text">The MIBiG ID of the BGC containing the enzyme</div>
        </div>
    `;
}

function createHtmlRef(data, className, jsonID) {
    return `
        <div class="${className}">
            <div class="row d-flex align-items-center mb-1">
                <div class="col">
                    <div class="form-floating">
                        <input type="text" name="${jsonID}" id="${jsonID}" class="form-control" aria-describedby="${jsonID}Help"  value="${data}" required>
                        <label for="${jsonID}" class="form-label">Reference DOI</label>
                        <div id="${jsonID}Help" class="form-text">A Digital Object Identifier (DOI), should begin with 'doi:'</div>
                    </div>
                </div>
                <div class="col-auto mx-auto">
                    <button type="button" class="btn btn-danger" onclick="removeField(this, '.${className}')">Remove</button>
                </div>
            </div>
        </div>
    `;
}

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




function createHtmlReaction(data = {}, index) {
    const reactionEntry = document.createElement('div')
    reactionEntry.classList.add('reaction')
    reactionEntry.innerHTML = `
        <div class="reaction">
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
                                    <div id="reaction[${index}]smarts-field"></div>
                                </div>
                            </div>
                            <div class="row mb-2">
                                <div class="col">
                                    <div id="reaction[${index}]description-field"></div>
                                </div>
                            </div>
                            <div class="row mb-2">
                                <div class="col">
                                    <div id="reaction[${index}]rhea-field"></div>
                                </div>
                                <div class="col">
                                    <div id="reaction[${index}]ec-field"></div>
                                </div>
                            </div>
                            <div class="row mb-2">
                                <div class="col">
                                    <div id='"reaction[${index}]tailoring-field"></div>
                                </div>
                            </div>
                            <div class="row mb-2">
                                <div class="col">
                                    <div id="reaction[${index}]evidencecode-field"></div>
                                </div>
                            </div>
                            <div class="row mb-2">
                                <div class="col">
                                    <div id="reaction[${index}]ref-field"></div>
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
                        </div>
                    </div>
                </div>
            </div>
        </div>
    `;

    const container = document.getElementById('reactions');
    container.appendChild(reactionEntry);

    const smartsHtml = addHtmlSmarts(data, index);
    const containerSmarts = document.getElementById(`reaction[${index}]smarts-field`);
    containerSmarts.insertAdjacentHTML('beforeend', smartsHtml);

//    TODO(MMZ 5.11.) Continue here: append html to the places where they should be

}

function addHtmlSmarts(data, index) {
    return `
        <div class="form-floating">
            <input id="reaction[${index}]smarts" name="reaction[${index}]smarts" class="form-control" aria-describedby="ReactionSmartsHelp" type="text" value='${data?.reactionSMARTS ?? ""}' required>
            <label for="enzyme_name" class="form-label">Reaction SMARTS</label>
            <div id="ReactionSmartsHelp" class="form-text">The reaction SMARTS depicting the reaction</div>
        </div>
    `;
}
