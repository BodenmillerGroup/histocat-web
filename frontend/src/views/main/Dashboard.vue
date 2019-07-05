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
  import { Component, Vue } from 'vue-property-decorator';
  import { readExperiments, readTags } from '@/modules/experiment/getters';
  import { dispatchGetExperiments, dispatchGetTags } from '@/modules/experiment/actions';
  import ExperimentCard from '@/components/ExperimentCard.vue';

  @Component({
    components: { ExperimentCard },
  })
  export default class Dashboard extends Vue {
    tags: string[] = [];

    get items(): any[] {
      const list = readTags(this.$store);
      return list.map((item) => {
        return {
          text: item,
        };
      });
    }

    get experiments() {
      const all = readExperiments(this.$store);
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
      await dispatchGetExperiments(this.$store);
      await dispatchGetTags(this.$store);
    }
  }
</script>
