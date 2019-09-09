<template>
  <v-banner
    v-if="!activeDataset || !activeAcquisition"
    icon="mdi-alert-circle-outline"
  >
    Please select acquisition and dataset
  </v-banner>
  <v-row
    v-else
    no-gutters
    class="chart-container"
  >
    <v-col
      :cols="columns"
    >
      <v-chart
        :options="options"
        autoresize
      />
    </v-col>
    <v-col
      v-if="showOptions"
      cols="3"
    >
      <v-card tile>
        <v-card-title>t-SNE Settings</v-card-title>
        <v-card-text>
          <v-chip-group
            v-model="selectedItems"
            multiple
            column
            active-class="primary--text"
          >
            <v-chip
              v-for="item in items"
              :key="item"
              :value="item"
              small
            >
              {{ item }}
            </v-chip>
          </v-chip-group>
          <v-card-actions class="input-row">
            <v-btn
              @click="selectAll"
              small
              :disabled="selectedItems.length === items.length"
            >
              Select all
            </v-btn>
            <v-btn
              @click="clearAll"
              small
              :disabled="selectedItems.length === 0"
            >
              Clear all
            </v-btn>
          </v-card-actions>
          <v-select
            class="input-row"
            :items="componentsNumbers"
            v-model.number="componentNumber"
            label="Number of components"
            hide-details
          ></v-select>
          <v-select
            class="input-row"
            :items="heatmaps"
            v-model="heatmap"
            label="Heatmap"
            hint="Heatmap marker"
            item-text="label"
            item-value="value"
            persistent-hint
            clearable
          ></v-select>
        </v-card-text>
        <v-card-actions>
          <v-btn
            @click="submit"
            color="primary"
            block
            :disabled="selectedItems.length === 0"
          >
            Analyze
          </v-btn>
        </v-card-actions>
        <v-card-text>
          <v-select
            class="input-row"
            :items="results"
            v-model="result"
            label="Results"
            hint="t-SNE processed data"
            persistent-hint
            clearable
          ></v-select>
        </v-card-text>
        <v-card-actions>
          <v-btn
            @click="display"
            color="primary"
            block
            :disabled="!result"
          >
            Display
          </v-btn>
        </v-card-actions>
      </v-card>
    </v-col>
  </v-row>
</template>

<script lang="ts">
  import { analysisModule } from '@/modules/analysis';
  import { ITSNEData } from '@/modules/analysis/models';
  import { datasetModule } from '@/modules/datasets';
  import { experimentModule } from '@/modules/experiment';
  import { mainModule } from '@/modules/main';
  import { settingsModule } from '@/modules/settings';
  import { required } from '@/utils';
  import * as echarts from 'echarts';
  import 'echarts-gl';
  import 'echarts/lib/chart/line';
  import 'echarts/lib/chart/scatter';
  import 'echarts/lib/component/toolbox';
  import 'echarts/lib/component/tooltip';
  import 'echarts/lib/component/visualMap';
  import { Component, Vue, Watch } from 'vue-property-decorator';

  const commonOptions: echarts.EChartOption = {
    title: {
      text: 't-Distributed Stochastic Neighbor Embedding',
      left: 'center',
      top: 0,
    },
    animation: false,
    tooltip: {
      show: true,
    },
    toolbox: {
      show: true,
      right: '9%',
      feature: {
        restore: {
          show: true,
          title: 'Reset',
        },
        saveAsImage: {
          show: true,
          title: 'Export',
        },
      },
    },
  };

  @Component
  export default class TSNETab extends Vue {
    readonly mainContext = mainModule.context(this.$store);
    readonly experimentContext = experimentModule.context(this.$store);
    readonly datasetContext = datasetModule.context(this.$store);
    readonly analysisContext = analysisModule.context(this.$store);
    readonly settingsContext = settingsModule.context(this.$store);

    readonly required = required;
    readonly componentsNumbers: number[] = [2, 3];

    options: echarts.EChartOption = {};

    selectedItems: any[] = [];
    componentNumber = 2;
    heatmap: string | null = null;

    result: string | null = null;

    get heatmaps() {
      return this.activeDataset && this.activeDataset.input['neighbors_columns'] ? this.activeDataset.input['neighbors_columns'].map(item => item.substring(10, item.length)) : [];
    }

    get showOptions() {
      return this.mainContext.getters.showOptions;
    }

    get columns() {
      return this.showOptions ? 9 : 12;
    }

    get activeAcquisition() {
      return this.experimentContext.getters.activeAcquisition;
    }

    get activeDataset() {
      return this.datasetContext.getters.activeDataset;
    }

    get items() {
      return this.activeDataset && this.activeDataset.input['channel_map'] ? Object.keys(this.activeDataset.input['channel_map']) : [];
    }

    get results() {
      return this.activeDataset && this.activeDataset.output && this.activeDataset.output['tsne'] ? Object.keys(this.activeDataset.output['tsne']).map(item => item) : [];
    }

    selectAll() {
      this.selectedItems = this.items;
    }

    clearAll() {
      this.selectedItems = [];
    }

    async submit() {
      if (await this.$validator.validateAll()) {
        if (!this.activeDataset) {
          self.alert('Please select a dataset');
          return;
        }

        if (!this.activeAcquisition) {
          self.alert('Please select an acquisition');
          return;
        }

        await this.analysisContext.actions.submitTSNE({
          dataset_id: this.activeDataset.id,
          acquisition_id: this.activeAcquisition.id,
          n_components: this.componentNumber,
          markers: this.selectedItems,
          heatmap: this.heatmap ? `Neighbors_${this.heatmap}` : '',
        });
      }
    }

    async display() {
      if (!this.activeDataset) {
        self.alert('Please select a dataset');
        return;
      }

      if (!this.result) {
        self.alert('Please select result data');
        return;
      }

      await this.analysisContext.actions.getTSNEResult({
        datasetId: this.activeDataset.id,
        name: this.result ? this.result : '',
      });
    }

    get tsneData() {
      return this.analysisContext.getters.tsneData;
    }

    @Watch('tsneData')
    tsneDataChanged(data: ITSNEData) {
      if (data) {
        if (data.z) {
          this.plot3D(data);
        } else {
          this.plot2D(data);
        }
      }
    }

    private plot2D(data: ITSNEData) {
      const points = data.heatmap ?
        data.x.data.map((x, i) => {
          return [x, data.y.data[i], data.heatmap!.data[i]];
        }) :
        data.x.data.map((x, i) => {
          return [x, data.y.data[i]];
        });

      const options: echarts.EChartOption = {
        ...commonOptions,
        xAxis: {
          type: 'value',
          name: data.x.label,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        yAxis: {
          type: 'value',
          name: data.y.label,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        dataset: {
          source: points,
          dimensions: [
            { name: data.x.label, type: 'float' },
            { name: data.y.label, type: 'float' },
          ],
        },
        series: [
          {
            type: 'scatter',
            name: 'Scatter2D',
            symbolSize: 4,
            large: !data.heatmap,
            encode: {
              x: data.x.label,
              y: data.y.label,
              tooltip: [data.x.label, data.y.label],
            },
          },
        ],
      };

      if (data.heatmap) {
        (options.dataset as any).dimensions.push(
          { name: data.heatmap.label, type: 'float' },
        );

        options.visualMap = this.getVisualMap(data);
      }

      this.options = options;
    }

    private plot3D(data: ITSNEData) {
      const options = {
        ...commonOptions,
        grid3D: {},
        xAxis3D: {
          type: 'value',
          name: data.x.label,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        yAxis3D: {
          type: 'value',
          name: data.y.label,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        zAxis3D: {
          type: 'value',
          name: data.z!.label,
          nameTextStyle: {
            fontWeight: 'bold',
          },
        },
        dataset: {
          source: [
            data.x.data,
            data.y.data,
            data.z!.data,
          ],
          dimensions: [
            { name: data.x.label, type: 'float' },
            { name: data.y.label, type: 'float' },
            { name: data.z!.label, type: 'float' },
          ],
        },
        series: [
          {
            type: 'scatter3D',
            name: 'Scatter3D',
            seriesLayoutBy: 'row',
            symbolSize: 4,
            encode: {
              x: data.x.label,
              y: data.y.label,
              z: data.z!.label,
              tooltip: [data.x.label, data.y.label, data.z!.label],
            },
          },
        ],
      } as echarts.EChartOption;

      if (data.heatmap) {
        (options.dataset as any).dimensions.push(
          { name: data.heatmap.label, type: 'float' },
        );

        (options.dataset as any).source.push(
          data.heatmap.data,
        );

        options.visualMap = this.getVisualMap(data);
      }

      this.options = options;
    }

    private getVisualMap(data: ITSNEData): echarts.EChartOption.VisualMap[] {
      const min = Math.min(...data.heatmap!.data);
      const max = Math.max(...data.heatmap!.data);
      return [
        {
          type: 'continuous',
          orient: 'horizontal',
          left: 'center',
          text: ['Max', 'Min'],
          calculable: true,
          realtime: false,
          min: min,
          max: max,
          inRange: {
            color: ['#4457cc', '#f45c00'],
          },
        },
      ];
    }
  }
</script>

<style scoped>
  .chart-container {
    height: calc(100vh - 154px);
  }

  .input-row {
    margin-bottom: 32px;
  }
</style>
