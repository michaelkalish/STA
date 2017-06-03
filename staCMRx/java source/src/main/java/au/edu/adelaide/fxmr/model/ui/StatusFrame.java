package au.edu.adelaide.fxmr.model.ui;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.concurrent.atomic.AtomicBoolean;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;
import java.awt.Font;
import au.edu.adelaide.fxmr.model.SolverListener;

/**
 * A Swing implementation of a cancellable status dialog. This is mainly in
 * place to allow tasks to be cancelled while they are running in MATLAB.
 * 
 * 
 */
public class StatusFrame extends JFrame implements ActionListener, SolverListener {
	private static final long serialVersionUID = 1L;
	private AtomicBoolean running;
	private JButton cancelButton = new JButton("Cancel");
	private JTextArea statusTextArea = new JTextArea();
	private final StringBuilder sb = new StringBuilder();

	public StatusFrame() {
		this("jCMRx");
	}

	public StatusFrame(String title) {
		this(title, new AtomicBoolean(true));
	}
	
	public StatusFrame(String title, AtomicBoolean running) {
		super(title + " " + version());
		this.running = running;
		cancelButton.addActionListener(this);
		statusTextArea.setEnabled(false);
		statusTextArea.setDisabledTextColor(Color.BLACK);
		statusTextArea.setFont(new Font("Courier New", Font.PLAIN, 12));

		JPanel panel = new JPanel(new BorderLayout(5, 5));

		panel.setBorder(BorderFactory.createEmptyBorder(5, 5, 5, 5));

		setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);

		setSize(540, 420);
		panel.add(new JScrollPane(statusTextArea), BorderLayout.CENTER);
		panel.add(cancelButton, BorderLayout.SOUTH);

		add(panel);

		setVisible(true);
	}

	public static String version() {
		String version = StatusFrame.class.getPackage().getImplementationVersion();
		if (version == null)
			return "dev";
		return version;
	}

	@Override
	public void actionPerformed(ActionEvent e) {
		running.set(false);
		statusTextArea.setText("Cancelling...");
	}

	public boolean updateStatus(final String status) {
		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				statusTextArea.setText(status);
			}
		});
		return running.get();
	}

	@Override
	public void setFinished() {
		setVisible(false);
	}

	@Override
	public boolean updateStatus(double fFloor, double fBar, double upperFloor, int size, int[] nIterThread, int collisions, int fBarReductions,
			int cyclesAvoided) {
		sb.setLength(0);
		sb.append("fFloor =\t");
		sb.append(fFloor);
		sb.append("\t-\t");
		sb.append(upperFloor);
		sb.append("\nfBar =\t");
		sb.append(fBar);
		sb.append("\nL size =\t");
		sb.append(size);
		sb.append("\n");

		int sum = 0;
		int nThread = nIterThread.length;
		for (int i = 0; i < nThread; i++) {
			sb.append("\nThread ");
			sb.append(i);
			sb.append(" Iter =\t");
			sb.append(nIterThread[i]);
			sum += nIterThread[i];
		}
		sb.append("\n\nTotal Iter =\t");
		sb.append(sum);

		sb.append("\n\nCollisions=\t");
		sb.append(collisions);
		sb.append("\nfBar Reductions= ");
		sb.append(fBarReductions);
		sb.append("\nAvoided Cycles=\t");
		sb.append(cyclesAvoided);

		sb.append("\n\nMem = \t");
		sb.append((Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory()) / (1024.0 * 1024.0));
		sb.append("\nMax = \t");
		sb.append(Runtime.getRuntime().totalMemory() / (1024.0 * 1024.0));

		SwingUtilities.invokeLater(new Runnable() {
			public void run() {
				statusTextArea.setText(sb.toString());
			}
		});

		return running.get();
	}

	public boolean isRunning() {
		return running.get();
	}

	public void setCancelText(String text) {
		cancelButton.setText(text);
	}
}
