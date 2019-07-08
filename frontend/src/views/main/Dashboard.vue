<template>
  <v-flex>
    <v-select
      v-model="tags"
      :items="items"
      chips
      deletable-chips
      clearable
      label="Tags"
      multiple
      solo
      class="mt-3 mx-3"
    ></v-select>

    <masonry
      :cols="3"
      :gutter="30"
    >
      <ExperimentCard
        v-for="experiment in experiments"
        :experiment="experiment"
        :key="experiment.id"
      />
    </masonry>
  </v-flex>
</template>

<script lang="ts">
  import ExperimentCard from '@/components/ExperimentCard.vue';
  import { experimentModule } from '@/modules/experiment';
  import { Component, Vue } from 'vue-property-decorator';

  @Component({
    components: { ExperimentCard },
  })
  export default class Dashboard extends Vue {
    experimentContext = experimentModule.context(this.$store);

    tags: string[] = [];

    get items(): any[] {
      const list = this.experimentContext.getters.tags;
      return list.map((item) => {
        return {
          text: item,
        };
      });
    }

    get experiments() {
      const all = this.experimentContext.getters.experiments;
      if (this.tags.length === 0) {
        return all;
      } else {
        return all.filter((experiment) => {
          if (experiment.tags && this.tags.some(r => experiment.tags.includes(r))) {
            return experiment;
          }
        });
      }
    }

    async mounted() {
      await this.experimentContext.actions.getExperiments();
      await this.experimentContext.actions.getTags();
    }
  }
</script>
