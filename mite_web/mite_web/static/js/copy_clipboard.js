var clipboard = new ClipboardJS('.btn-light');

clipboard.on('success', function(e) {
     alert('Copied to clipboard!');
     e.clearSelection();
});

clipboard.on('error', function(e) {
     alert('Failed to copy text.');
});