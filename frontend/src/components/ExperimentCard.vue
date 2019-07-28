<template>
  <v-card
    tile
    class="ma-6 pa-1"
  >
    <v-card-title>
      <v-layout column>
        <h5 class="headline">{{experiment.name}}</h5>
        <span class="caption"><v-icon small>mdi-calendar-outline</v-icon> {{createdAt}}</span>
      </v-layout>
    </v-card-title>
    <v-card-text>
      {{experiment.description}}
    </v-card-text>
    <v-card-text>
      <v-chip
        :key="item"
        v-for="item in experiment.tags"
        label
        small
        class="mr-1"
      >
        <v-icon small left>mdi-tag-outline</v-icon>
        {{ item }}
      </v-chip>
    </v-card-text>
    <v-card-actions>
      <v-btn tile color="indigo" dark :to="{name: 'main-experiment', params: {id: experiment.id}}">Open</v-btn>
      <v-spacer></v-spacer>
      <v-tooltip bottom>
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on" :to="{name: 'main-experiment-share', params: {id: experiment.id}}">
            <v-icon>mdi-share-variant</v-icon>
          </v-btn>
        </template>
        <span>Share experiment</span>
      </v-tooltip>
    </v-card-actions>
  </v-card>
</template>

<script lang="ts">
  import { Component, Prop, Vue } from 'vue-property-decorator';
  import { IExperiment } from '@/modules/experiment/models';

  @Component
  export default class ExperimentCard extends Vue {
    @Prop(Object) experiment?: IExperiment;

    get createdAt() {
      return this.experiment && new Date(this.experiment.created_at).toUTCString();
    }
  }
</script>
