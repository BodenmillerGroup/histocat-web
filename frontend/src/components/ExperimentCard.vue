<template>
  <v-card class="ma-3 pa-3" tile hover max-width="400">
    <v-card-title>
      <v-layout column>
        <h5 class="headline">{{experiment.name}}</h5>
        <span class="caption">{{createdAt}}</span>
      </v-layout>
    </v-card-title>
    <v-card-text>
      {{experiment.description}}
    </v-card-text>
    <v-card-actions>
      <v-btn :to="{name: 'main-experiment', params: {id: experiment.id}}">View</v-btn>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';
  import { IExperiment } from '@/modules/experiment/models';
  import { DateTime } from 'luxon';

  @Component
  export default class ExperimentCard extends Vue {
    @Prop(Object) experiment?: IExperiment;

    get createdAt() {
      return this.experiment && DateTime.fromISO(this.experiment.created_at).toFormat('yyyy-MM-dd HH:mm:ss');
    }
  }
</script>
