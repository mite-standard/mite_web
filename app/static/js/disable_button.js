 // prevents repeated POST requests
function disableButton(form, button) {
  form.addEventListener('submit', function (event) {
      if (!form.checkValidity()) {
          return;
      }
      setTimeout(() => {
            button.disabled = true;
            button.innerText = 'Submitting...';
        }, 0);
  });
}
