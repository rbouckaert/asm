package asm.tools;

import java.io.File;
import java.io.IOException;
import java.net.URL;

import beast.base.core.ProgramStatus;
import beastfx.app.inputeditor.MyAction;
import beastfx.app.tools.LogAnalyser;
import beastfx.app.util.Alert;
import beastfx.app.util.FXUtils;
import javafx.application.Application;
import javafx.beans.property.SimpleIntegerProperty;
import javafx.beans.property.SimpleStringProperty;
import javafx.event.ActionEvent;
import javafx.scene.Cursor;
import javafx.scene.Scene;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Menu;
import javafx.scene.control.MenuBar;
import javafx.scene.control.MenuItem;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.scene.input.KeyCode;
import javafx.scene.input.KeyCodeCombination;
import javafx.scene.input.KeyCombination;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;

public class Monitor extends Application {
	LineChart<Number, Number> chart;
	static LineChart.Series<Number, Number> series;
	static LineChart.Series<Number, Number> series2;
	static int traceIndex = 1;
	static String label, filename;
	static double lower, upper;
	static LogAnalyser trace;
	TableView<LabelInfo> table;

	public class LabelInfo {

		public SimpleStringProperty name;
		public SimpleIntegerProperty burnin;
		private int index;

		LabelInfo(String name, int burnin, int index) {
			this.index = index;
			this.name = new SimpleStringProperty(name);
			this.burnin = new SimpleIntegerProperty(burnin);
		}

		public String getLabel() {
			return name.get();
		}

		public void setLabel(String fname) {
			name.set(fname);
		}

		public void setBurnin(Integer burnin) {
			this.burnin.set(burnin);
		}

		public Integer getBurnin() {
			return burnin.get();
		}

		public SimpleStringProperty getNameProperty() {
			return name;
		}

		public SimpleIntegerProperty getBIntegerProperty() {
			return burnin;
		}
		
		public int getIndex() {
			return index;
		}

	};

	@Override
	public void start(final Stage stage) throws Exception {
		NumberAxis xAxis = new NumberAxis();
		xAxis.setForceZeroInRange(false);
		xAxis.setLabel(filename);
		NumberAxis yAxis = new NumberAxis();
		yAxis.setForceZeroInRange(false);
		yAxis.setLabel(label);
		yAxis.setUpperBound(upper);
		yAxis.setLowerBound(lower);
		yAxis.setAutoRanging(false);
		chart = new LineChart<>(xAxis, yAxis);
		chart.setLegendVisible(false);
		chart.setCreateSymbols(false);
		chart.getData().add(series);
		chart.getData().add(series2);

		MenuBar menuBar = createMenuBar();

		BorderPane pane = new BorderPane(chart);
		pane.setTop(menuBar);
		table = new TableView<>();
		TableColumn<LabelInfo, String> col = new TableColumn<>("Label");
		col.setCellValueFactory(new PropertyValueFactory<LabelInfo, String>("Label"));
		col.setPrefWidth(100);
		table.getColumns().add(col);
		TableColumn<LabelInfo, String> col2 = new TableColumn<>("Burnin");
		col2.setPrefWidth(50);
		col2.setCellValueFactory(new PropertyValueFactory<LabelInfo, String>("Burnin"));
		table.getColumns().add(col2);
		initTable();
		pane.setLeft(table);

		Scene scene = new Scene(pane, 1024, 600);
		URL url = Monitor.class.getResource("style.css");
		String css = url.toExternalForm();
		scene.getStylesheets().add(css);
		stage.setScene(scene);
		stage.show();
	}

	class ActionLoad extends MyAction {

		public ActionLoad() {
			super("Load", "Load Trace Log File", "open",
					new KeyCodeCombination(KeyCode.O, KeyCombination.SHORTCUT_DOWN));
		}

		@Override
		public void actionPerformed(ActionEvent ae) {
			File file = FXUtils.getLoadFile("Load Trace Log File", new File(ProgramStatus.g_sDir), "Trace log files",
					"log");
			if (file != null) {
				chart.getScene().getRoot().setCursor(Cursor.WAIT);
				try {
					trace = new LogAnalyser(file.getAbsolutePath());
					initTable();
					updateChart();
					chart.getXAxis().setLabel(file.getName());
					chart.getScene().getRoot().setCursor(Cursor.DEFAULT);
				} catch (IOException e) {
					e.printStackTrace();
					Alert.showMessageDialog(null, "Something went wrong loading the file: " + e.getMessage());
				}
			}
		}
	} // ActionLoad

	class ActionQuit extends MyAction {
		public ActionQuit() {
			super("Exit", "Exit Program", "exit", new KeyCodeCombination(KeyCode.Q, KeyCombination.SHORTCUT_DOWN));
		} // c'tor

		@Override
		public void actionPerformed(ActionEvent ae) {
			System.exit(0);
		}
	} // class ActionQuit

	private MenuBar createMenuBar() {
		MenuBar menuBar = new MenuBar();
		Menu fileMenu = new Menu("_File");
		fileMenu.setId("File");
		fileMenu.setMnemonicParsing(true);
		menuBar.getMenus().add(fileMenu);
		fileMenu.getItems().clear();

		MenuItem a_load = new ActionLoad();
		MenuItem a_quit = new ActionQuit();

		fileMenu.getItems().add(a_load);
		fileMenu.getItems().add(a_quit);
		return menuBar;
	}

	private void initTable() {
		table.getItems().clear();
		for (int i = 0; i < trace.getLabels().size(); i++) {
			Double[] t = trace.getTrace(i + 1);
			LabelInfo e = new LabelInfo(trace.getLabels().get(i), burnIn(t, t.length), i + 1);
			table.getItems().add(e);
		}

		table.getSelectionModel().selectedItemProperty().addListener((observable, oldValue, newValue) -> {
			updateChart();
		});
		table.getSelectionModel().select(0);

	}

	private void updateChart() {

		int selected = table.getSelectionModel().getSelectedIndex();
		if (selected >= trace.getLabels().size()) {
			selected = 0;
		}
		Double[] t = trace.getTrace(selected + 1);
		// series = new LineChart.Series<>();
		// series.getData().clear();
		for (int i = 0; i < t.length; i++) {
			if (i < series.getData().size()) {
				series.getData().set(i, new XYChart.Data<Number, Number>(i, t[i]));
			} else {
				series.getData().add(new XYChart.Data<Number, Number>(i, t[i]));
			}
		}
		while (series.getData().size() > t.length) {
			series.getData().remove(series.getData().size() - 1);
		}

		int burnin = burnIn(t, t.length);
		// series2 = new LineChart.Series<>();
		lower = t[burnin];
		upper = lower;
		for (int i = burnin; i < t.length; i++) {
			if (i - burnin < series2.getData().size()) {
				series2.getData().set(i - burnin, new XYChart.Data<Number, Number>(i, t[i]));
			} else {
				series2.getData().add(new XYChart.Data<Number, Number>(i, t[i]));
			}
			lower = Math.min(lower, t[i]);
			upper = Math.max(upper, t[i]);
		}
		while (series2.getData().size() > t.length - burnin) {
			series2.getData().remove(series2.getData().size() - 1);
		}
		// chart.getData().clear();
		// chart.getData().set(0,series);
		// chart.getData().set(1,series2);

		try {
			label = trace.getLabels().get(selected);
			NumberAxis yAxis = (NumberAxis) chart.getYAxis();
			yAxis.setLabel(label);
			yAxis.setUpperBound(upper);
			yAxis.setLowerBound(lower);
			yAxis.setAutoRanging(false);
		} catch (IndexOutOfBoundsException e) {
			// ignore
		}
	}

	private static int burnIn(Double[] trace, int end) {
		int WINDOW_SIZE = 10;
		// calc mean and stdev of last 25%
		double m2 = 0, sq2 = 0;
		int lb = 3 * end / 4;
		for (int i = lb; i < end; i++) {
			double d = trace[i];
			m2 += d;
		}
		m2 = m2 / (end - lb);
		for (int i = lb; i < end; i++) {
			double d = trace[i];
			sq2 += (d - m2) * (d - m2);
		}
		double stdev2 = Math.sqrt(sq2 / (end - lb - 1));

		// pick first moving average that fits inside the m2 +/- stdev2 range
		for (int i = 0; i + WINDOW_SIZE < lb; i++) {
			double sum = 0;
			for (int j = 0; j < WINDOW_SIZE; j++) {
				sum += trace[i + j];
			}
			sum /= WINDOW_SIZE;
			if (m2 - stdev2 < sum && sum < m2 + stdev2) {
				return i + WINDOW_SIZE;
			}
		}
		return lb;
	}

	public static void main(String[] args) throws IOException {
		filename = args[0];
		trace = new LogAnalyser(filename, 0, true, false);

		if (args.length > 1) {
			for (int i = 0; i < trace.getLabels().size(); i++) {
				if (args[1].equals(trace.getLabels().get(i))) {
					traceIndex = i + 1;
				}
			}
		}
		label = trace.getLabels().get(traceIndex - 1);
		Double[] t = trace.getTrace(traceIndex);
		series = new LineChart.Series<>();
		for (int i = 0; i < t.length; i++) {
			series.getData().add(new XYChart.Data<Number, Number>(i, t[i]));
		}

		int burnin = burnIn(t, t.length);
		series2 = new LineChart.Series<>();
		lower = t[burnin];
		upper = lower;
		for (int i = burnin; i < t.length; i++) {
			series2.getData().add(new XYChart.Data<Number, Number>(i, t[i]));
			lower = Math.min(lower, t[i]);
			upper = Math.max(upper, t[i]);
		}

		launch(Monitor.class, args);
	}

}